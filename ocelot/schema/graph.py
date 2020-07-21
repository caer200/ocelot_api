import itertools
import warnings
import hashlib
from collections import OrderedDict
from operator import eq
from ocelot.routines.fileop import stringkey, intkey
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from networkx import Graph
from pymatgen.core.periodic_table import Element
from pymatgen.vis.structure_vtk import EL_COLORS
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from ocelot.routines.conformerparser import ACParser
from ocelot.schema.rdfunc import RdFunc

"""
a graph of integer nodes with a "symbol" attribute
all edges are defined by nodes, no attribute

molgraph is a molecule with single bonds only
only node-induced subgraph in substructure search

If G’=(N’,E’) is a node-induced subgraph, then:
N’ is a subset of N E’ is the subset of edges in E relating nodes in N’
"""


def to_fracrgb(t):
    return [x / 255 for x in t]

class GraphError(Exception): pass

def hashstr2int(s, hashmethod = hashlib.sha256):
    return int(hashmethod(s.encode('utf-8')).hexdigest(), 16) % 10**8


class BasicGraph:

    def __init__(self, graph: nx.Graph):
        self.graph = graph
        # self.rings = nx.minimum_cycle_basis(self.graph)  # technically sssr
        # self.rings = sorted(self.rings, key=lambda x: len(x))
        # self.nrings = len(self.rings)

    def __len__(self):
        return len(self.graph.nodes)

    # TODO find a better way to hash a graph schema
    def __hash__(self):
        return self.hash_nxgraph(self.graph)
        # return self.hash_via_smiles(self)

    @property
    def props_for_hash(self):
        """
        see https://stackoverflow.com/questions/46999771/
        use with caution...
        """
        t = nx.triangles(self.graph)
        c = nx.number_of_cliques(self.graph)
        ele = nx.get_node_attributes(self.graph, 'symbol')
        dv = self.graph.degree
        props = [(dv[v], t[v], c[v], ele[v]) for v in self.graph]
        props.sort()
        return props

    @staticmethod
    def hash_via_smiles(mg):
        smiles, rdmol, _, _ = mg.to_rdmol()
        # return hash("{}: {}".format(mg.__class__.__name__, smiles))
        return hashstr2int(smiles)

    @staticmethod
    def hash_nxgraph(g: nx.Graph):
        t = nx.triangles(g)
        c = nx.number_of_cliques(g)
        ele = nx.get_node_attributes(g, 'symbol')
        dv = g.degree
        props = [(dv[v], t[v], c[v], ele[v]) for v in g]
        props.sort()
        s = str(props)
        return hashstr2int(s)
        # return hash(tuple(props))
        # pnGraph=pn.Graph(g.number_of_nodes());
        # edg=list(g.edges);
        # for E in edg:
        #     pnGraph.connect_vertex(E[0],E[1]);
        # return hash(pn.certificate(pnGraph));

    def __repr__(self):
        outs = ["{}:".format(self.__class__.__name__)]
        for n in self.graph.nodes:
            outs.append('{} {}'.format(n, self.symbols[n]))
        return "; ".join(outs)

    def __eq__(self, other):
        """
        equality is measured based on graph isomorphism if smiles hash failed
        """
        try:
            return hash(self) == hash(other)
        except:
            g1 = self.graph
            g2 = other.graph
            return nx.is_isomorphic(g1, g2, node_match=iso.generic_node_match('symbol', None, eq))

    def graph_similarity(self, other):
        """
        calculate graph similarity (GED) between two MolGraphs

        :param other:
        :return:
        """
        return nx.algorithms.graph_edit_distance(self.graph, other.graph,
                                                 node_match=iso.generic_node_match('symbol', None, eq))

    def is_subgraph(self, other):
        """
        True if other has a node-induced subgraph is isomorphic to self
        """
        gm = iso.GraphMatcher(other.graph, self.graph, node_match=iso.generic_node_match('symbol', None, eq))
        return gm.subgraph_is_isomorphic()

    @property
    def symbols(self):
        """
        a dict s.t. symbols[node] gives symbol
        """
        return nx.get_node_attributes(self.graph, 'symbol')

    @classmethod
    def from_smiles(cls, smiles:str):
        m = Chem.MolFromSmiles(smiles)
        if 'H' in smiles or 'h' in smiles:
            hm = m
        else:
            hm = Chem.AddHs(m)
        return cls.from_rdmol(hm)

    @classmethod
    def from_rdmol(cls, rdmol: Chem.Mol, atomidx2nodename=None):
        """
        :param rdmol:
        :param atomidx2nodename: the dict to convert atomidx in rdmol to graph node, atomidx2nodename[rdmolid] == nodename
        if not None then nodename will be set just based on atomidx
        :return:
        """
        g = nx.Graph()
        if atomidx2nodename:
            index_dict_no = atomidx2nodename
            for atom in rdmol.GetAtoms():
                g.add_node(index_dict_no[atom.GetIdx()], symbol=atom.GetSymbol(), )
            for bond in rdmol.GetBonds():
                g.add_edge(
                    index_dict_no[bond.GetBeginAtomIdx()],
                    index_dict_no[bond.GetEndAtomIdx()],
                )
        else:
            for atom in rdmol.GetAtoms():
                g.add_node(atom.GetIdx(),
                           symbol=atom.GetSymbol(),
                           )
            for bond in rdmol.GetBonds():
                g.add_edge(
                    bond.GetBeginAtomIdx(),
                    bond.GetEndAtomIdx(),
                )
        if not nx.is_connected(g):
            raise GraphError('the graph is not connected!')
        return cls(g)

    @classmethod
    def from_dict(cls, d: dict):
        g = nx.from_dict_of_lists(intkey(d['graph_dict_of_lists']))
        nx.set_node_attributes(g, intkey(d['graph_node_symbols']), name='symbol')
        return cls(g)

    def to_rdmol(self, sani=True, charge=0, charged_fragments=None, force_single=False, expliciths=True):
        """
        convert to rdmol, it will first try to convert to a radical, if failed b/c of rdkit valence rules,
        it will try to convert to charged fragment

        :param sani: switch to call default sanitizer of rdkit
        :param int charge: molecule charge
        :param bool charged_fragments: switch to assign atomic charge
        :param bool force_single: switch to force all single bonds in the resulting molecule
        :param bool expliciths: add hydrogen explicitly
        :return: a rdmol object, canonical smiles, atomidx2nodename, nodename2atomidx
        """
        # I cannot use smilestomol as I need to know how atom ids are mapped
        nodes = sorted(list(self.graph.nodes(data=True)), key=lambda x: x[0])
        atomidx2nodename = {}  # d[atom idx in the new rdmol] = original graph node
        nodename2atomidx = {}


        atom_number_list = []
        for i in range(len(nodes)):
            graph_node = nodes[i][0]  # this could be different form i!
            symbol = nodes[i][1]['symbol']
            z = Element(symbol).Z
            atom_number_list.append(z)
            atomidx2nodename[i] = graph_node
            nodename2atomidx[graph_node] = i

        adj = nx.convert.to_dict_of_dicts(self.graph)

        new_ac = np.zeros((len(nodes), len(nodes))).astype(int)  # atomic connectivity
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if atomidx2nodename[j] in adj[atomidx2nodename[i]].keys():
                    new_ac[i, j] = 1
                    new_ac[j, i] = 1

        apriori_radicals = {}
        if hasattr(self, 'joints'):  # LBYL
            for k in self.joints.keys():
                apriori_radicals[nodename2atomidx[k]] = len(self.joints[k])
        else:
            apriori_radicals = None

        ap = ACParser(sani=sani, ac=new_ac, atomnumberlist=atom_number_list, charge=charge, apriori_radicals=apriori_radicals)
        if charged_fragments is None:
            try:
                mol, smiles = ap.parse(charged_fragments=False, force_single=force_single, expliciths=expliciths)
            except Chem.rdchem.AtomValenceException:
                warnings.warn('AP parser cannot use radical scheme, trying to use charged frag')
                mol, smiles = ap.parse(charged_fragments=True, force_single=force_single, expliciths=expliciths)
        else:
            mol, smiles = ap.parse(charged_fragments=charged_fragments, force_single=force_single, expliciths=expliciths)
        return mol, smiles, atomidx2nodename, nodename2atomidx

    def as_dict(self):
        d = OrderedDict()
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['graph_node_symbols'] = stringkey(self.symbols)
        d['graph_dict_of_lists'] = stringkey(nx.to_dict_of_lists(self.graph))
        # d['rings'] = self.rings
        # d['nrings'] = self.nrings
        return d

    def draw(self, output: str, dpi=600):
        """
        draw graph, using jmol color scheme for elements

        :param output:
        :param dpi:
        :return:
        """
        draw_chemgraph(self.graph, output, dpi)


def draw_chemgraph(graph: nx.Graph, output: str, dpi=600):
    f = plt.figure()
    symbols = nx.get_node_attributes(graph, 'symbol')
    colors = [to_fracrgb(EL_COLORS['Jmol'][symbols[n]]) for n in graph]
    nx.draw(graph, with_labels=True, labels=symbols, node_color=colors)
    f.savefig(output, dpi=dpi)


class FragmentGraph(BasicGraph):

    def __init__(self, graph: nx.Graph, joints: dict, partition_scheme: str):
        """
        fragment is a BasicGraph with joints

        I feel it's not necessary to add parent BasicGraph as an attribute, but I maybe wrong

        :param graph:
        :param joints: a dict, keys are nodes in the subgraph that have parent graph edges not in this subgraph
        """
        super().__init__(graph)
        self.partition_scheme = partition_scheme
        self.graph = graph
        self.joints = joints
        self.rings = nx.minimum_cycle_basis(self.graph)  # technically sssr
        self.rings = sorted(self.rings, key=lambda x: len(x))

        # we don't need ring graph in frag, that should be used in partition process
        # self.rings = nx.minimum_cycle_basis(self.graph)
        # self.rings = sorted(self.rings, key=lambda x: len(x))
        # self.nrings = len(self.rings)

    def as_dict(self):
        d = super().as_dict()
        d['joints'] = stringkey(self.joints)
        d['partition_scheme'] = self.partition_scheme
        return d

    @classmethod
    def from_dict(cls, d: dict):
        g = nx.from_dict_of_lists(intkey(d['graph_dict_of_lists']))
        nx.set_node_attributes(g, intkey(d['graph_node_symbols']), name='symbol')
        return cls(g, intkey(d['joints']), d['partition_scheme'])


class BackboneGraph(FragmentGraph):  # is this necessary?
    def __init__(self, graph, joints, partition_scheme='lgcr'):
        super().__init__(graph, joints, partition_scheme)
        # self.ringfp = set([len(r) for r in self.rings])

    # def as_dict(self):
    #     d = super().as_dict()
    #     # d['ringfp'] = self.ringfp
    #     return d


class SidechainGraph(FragmentGraph):
    def __init__(self, graph, joints, partition_scheme):
        """
        it is not possible for a side chain to have two frag-joints connected to one bone-joint, as this should be
        already in backbone

        side chain (frags) cannot have 2 side joints, each connected to one bone-joint as this would be a ring fused
         to backbone (side chain is a connected component)

        :var self.rankmap: this gives the rank for each node, rank is the shortest path length from the node to sc_joint
        """
        super().__init__(graph, joints, partition_scheme)
        if len(self.joints.keys()) > 1:
            # in general, sc from chrom partition should not be used
            # warnings.warn('sidechain has more than one side joint... this could happen when using chrom partition')
            raise GraphError('sidechain has more than one side joint... impossible!')
        if len(self.joints.keys()) == 0:
            raise GraphError('sidechain has no side joint, the molceule may not be a connected graph')
        sc_joint = list(self.joints.keys())[0]
        self.rankmap = {}
        for node in self.graph:
            rank = nx.shortest_path_length(self.graph, source=node, target=sc_joint)
            self.rankmap[node] = rank
            # self.graph.nodes[node]['rank'] = rank

        # if self.nrings > 0:
        #     self.hasring = True
        # else:
        #     self.hasring = False

    def as_dict(self):
        d = super().as_dict()
        # d['hasring'] = self.hasring
        d['rankmap'] = stringkey(self.rankmap)
        return d


class MolGraph(BasicGraph):

    def __init__(self, graph):
        """
        this is the schema for molecule with fused rings

        :param graph:
        """
        super().__init__(graph)
        self.rings = nx.minimum_cycle_basis(self.graph)  # technically sssr
        self.rings = sorted(self.rings, key=lambda x: len(x))
        self.nrings = len(self.rings)
        if self.nrings < 1:
            warnings.warn('you are init a MolGraph with no ring!')

        self.ring_graph = nx.Graph()
        for ij in itertools.combinations(range(self.nrings), 2):
            i, j = ij
            ri = self.rings[i]
            rj = self.rings[j]
            shared_nodes = set(ri) & set(rj)
            ni_connected_to_rj = [ni for ni in ri if len(set(self.graph.neighbors(ni)) & set(rj)) == 1 and ni not in rj]
            # we set len(shared_nodes) > 0 to avoid breaking spiro bicycle
            if len(shared_nodes) > 0 or len(ni_connected_to_rj) > 0:
                self.ring_graph.add_edge(i, j, nshare=len(shared_nodes), nconnect=len(ni_connected_to_rj))

        self.connected_rings = self.get_rings('nconnect', 1)  # for rubrene len(self.connected_rings[0]) == 8
        self.fused_rings = self.get_rings('nshare', 2)  # for rubrene len(self.fused_rings[0]) == 4
        # notice if two rings share2 then they must connect1
        try:
            self.lgcr = self.connected_rings[0]
        except IndexError:
            try:
                self.lgcr = [self.rings[0]]
            except IndexError:
                raise GraphError('no rings in the MolGraph!')

        try:
            self.lgfr = self.fused_rings[0]
        except IndexError:
            try:
                self.lgfr = [self.rings[0]]
            except IndexError:
                raise GraphError('no rings in the MolGraph!')

    @classmethod
    def from_basicgraph(cls, basicgraph: BasicGraph):
        return cls(basicgraph.graph)

    def get_rings(self, k='nconnect', threshold=1):
        """
        get connected rings between which # of edges is larger than a threshold

        :param k:
        :param threshold:
        :return:
        """
        edges = [(u, v) for (u, v, d) in self.ring_graph.edges(data=True) if d[k] >= threshold]
        subgraph = self.ring_graph.edge_subgraph(edges).copy()
        return [[self.rings[iring] for iring in c] for c in
                sorted(nx.connected_components(subgraph), key=len, reverse=True)]

    @staticmethod
    def get_joints_and_subgraph(subg1_nodes: [int], g: Graph):
        """
        we assume subg1 \cap subg2 == empty and subg1 + subg2 == g
        """
        sg1 = nx.Graph()
        sg2 = nx.Graph()  # this is the graph of anything other than bone nodes
        symbols = nx.get_node_attributes(g, 'symbol')
        for node in g:
            s = symbols[node]
            if node in subg1_nodes:
                sg1.add_node(node, symbol=s)
            else:
                sg2.add_node(node, symbol=s)

        joints_subg1_as_keys = {}
        joints_subg2_as_keys = {}
        for n, nbrs in g.adj.items():
            for nbr, d in nbrs.items():
                if n in subg1_nodes and nbr in subg1_nodes:
                    sg1.add_edge(n, nbr, **d)
                elif n not in subg1_nodes and nbr not in subg1_nodes:
                    sg2.add_edge(n, nbr, **d)
                elif n in subg1_nodes and nbr not in subg1_nodes:
                    if n not in joints_subg1_as_keys.keys():
                        joints_subg1_as_keys[n] = [nbr]
                    else:
                        joints_subg1_as_keys[n].append(nbr)
                elif n not in subg1_nodes and nbr in subg1_nodes:
                    if n not in joints_subg2_as_keys.keys():
                        joints_subg2_as_keys[n] = [nbr]
                    else:
                        joints_subg2_as_keys[n].append(nbr)

        sg1.graph['joints'] = joints_subg1_as_keys

        sg2_components = []
        for c in sorted(nx.connected_components(sg2), key=len, reverse=True):
            nodes = list(c)
            joints_in_frag = {}
            for j in joints_subg2_as_keys.keys():
                if j in nodes:
                    joints_in_frag[j] = joints_subg2_as_keys[j]
            frag_graph = nx.Graph(joints=joints_in_frag)  # joints_in_frag[frag_node] is a list of bone joints
            frag_graph.add_nodes_from((n, g.nodes[n]) for n in nodes)
            frag_graph.add_edges_from(
                (n, nbr, d) for n, nbrs in g.adj.items() if n in nodes for nbr, d in nbrs.items() if
                nbr in nodes)
            sg2_components.append(frag_graph)
        return joints_subg1_as_keys, joints_subg2_as_keys, sg1, sg2_components

    @staticmethod
    def get_bone_and_frags_from_nxgraph(bone_graph: nx.Graph, fragments: [nx.Graph], scheme="bm"):
        if scheme in ["lgcr", "lgfr"]:
            gb = BackboneGraph(bone_graph, bone_graph.graph['joints'], partition_scheme=scheme)
            scs = [SidechainGraph(frag_graph, frag_graph.graph['joints'], partition_scheme=scheme) for frag_graph in
                   fragments]
        else:
            gb = FragmentGraph(bone_graph, bone_graph.graph['joints'], partition_scheme=scheme)
            scs = [FragmentGraph(frag_graph, frag_graph.graph['joints'], partition_scheme=scheme) for frag_graph in
                   fragments]
        return gb, scs

    @staticmethod
    def get_chrom_and_frags_from_nxgraph(chrom_graph, fragments):
        cg = FragmentGraph(chrom_graph, chrom_graph.graph['joints'], 'chrom')
        fgs = [FragmentGraph(frag_graph, frag_graph.graph['joints'], 'chrom') for frag_graph in fragments]
        return cg, fgs

    def partition(self, bone_selection='bm', additional_ring_criteria=None, with_halogen=True):
        """
        parition the molecule into a backbone graph and a list of fragments (graphs)

        :param with_halogen:
        :param bone_selection:
        :param additional_ring_criteria:
            this should be a function to check whether a ring (a set of siteids) meets additional conditions
            # this does nothing for chrom scheme
        :return:
        """
        if bone_selection == 'lgfr':
            rings = self.lgfr
        elif bone_selection == 'lgcr':
            rings = self.lgcr
        elif bone_selection == 'chrom':
            try:
                mol, smiles, atomidx2nodename, nodename2atomidx = self.to_rdmol()
            except:
                raise GraphError('to_rdmol failed in chrom partition!')
            if with_halogen:
                cgs = RdFunc.get_conjugate_group_with_halogen(mol)
            else:
                cgs = RdFunc.get_conjugate_group(mol)
            try:
                chromol, aid_to_newid_chromol = cgs[0]
            except IndexError:
                raise GraphError('rdkit cannot find a chrom here!')
            aids_in_chromol = list(aid_to_newid_chromol.keys())
            nodes_in_chromol = [atomidx2nodename[aid] for aid in aids_in_chromol]
            nodes_not_in_chromol = [sid for sid in self.graph if sid not in nodes_in_chromol]
            chromol_joints, other_joints, chromolsg, sg_components = MolGraph.get_joints_and_subgraph(
                nodes_in_chromol, self.graph)
            return chromolsg, sg_components
        #
        elif bone_selection == 'bm':
            try:
                mol, smiles, atomidx2nodename, nodename2atomidx = self.to_rdmol()
            except:
                raise GraphError('to_rdmol failed in chrom partition!')
            
            bm_mol = MurckoScaffold.GetScaffoldForMol(mol)
            backbone_nodes = set(mol.GetSubstructMatch(bm_mol))
            backbone_nodes = set([atomidx2nodename[node] for node in backbone_nodes])
            joints_bone_as_keys, joints_frag_as_keys, bone_graph, fragments = MolGraph.get_joints_and_subgraph(
                backbone_nodes,
                self.graph)

            return bone_graph, fragments
            
        else:
            raise GraphError('selection scheme not implemented in graph partition: {}'.format(bone_selection))

        if additional_ring_criteria is None:
            bonerings = rings
        else:
            bonerings = []
            for r in rings:
                if additional_ring_criteria(r):
                    bonerings.append(r)

        backbone_nodes = set([item for sublist in bonerings for item in sublist])
        joints_bone_as_keys, joints_frag_as_keys, bone_graph, fragments = self.get_joints_and_subgraph(backbone_nodes,
                                                                                                       self.graph)

        return bone_graph, fragments

    def partition_to_bone_frags(self, bone_selection='bm', additional_criteria=None):
        """
        partition the molecule into backbone and a list of frags

        the list of frags is reversely sorted by graph size

        :param additional_criteria:
        :param bone_selection:
        :return:
        """
        bone_graph, fragments = self.partition(bone_selection, additional_criteria)
        gb, scs = self.get_bone_and_frags_from_nxgraph(bone_graph, fragments, scheme=bone_selection)
        # gb = BackboneGraph(bone_graph, bone_graph.graph['joints'], partition_scheme=bone_selection)
        # scs = [SidechainGraph(frag_graph, frag_graph.graph['joints'], partition_scheme=bone_selection) for frag_graph in
        #        fragments]
        # scs = sorted(scs, key=lambda x: len(x.graph), reverse=True)
        return gb, scs
