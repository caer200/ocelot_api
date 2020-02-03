from ocelot.schema.element import Element
from ocelot.routines.conformerparser import ACParser
import networkx as nx
from collections import OrderedDict
import itertools
import networkx.algorithms.isomorphism as iso
from rdkit.Chem import Descriptors
from rdkit import Chem
# import copy
from pymatgen.core.periodic_table import Element
import numpy as np
import warnings
from operator import eq
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

"""
we use graph as a shadow rdkit.mol class

functions should include
1. substructure search
2. similarity calculation
3. side chain definition/descriptors
4. loop finding
5. backbon definition/descriptors
6. from/to conformation
7. display
8. decomposition

#TODO match is quite slow
"""


class GraphMolecule:

    def __init__(self, graph, hsmiles):
        """
        this is the basic object for a connected molecule, only connectivity and symbols are stored

        a fused ring is considered as a set of cycles in which at least one pair share at least 2 nodes

        :param graph:
        :param hsmiles:
        """
        self.graph = graph
        self.hsmiles = hsmiles
        self.rings = nx.minimum_cycle_basis(self.graph)  # technically sssr
        self.rings = sorted(self.rings, key=lambda x: len(x))
        self.symbols = nx.get_node_attributes(self.graph, 'symbol')
        self.nrings = len(self.rings)

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
        if len(self.connected_rings) > 0:
            self.lgcr = self.connected_rings[0]
        else:
            if len(self.rings) > 0:
                self.lgcr = [self.rings[0]]
            else:
                self.lgcr = None
        if len(self.fused_rings) > 0:
            self.lgfr = self.fused_rings[0]
        else:
            if len(self.rings) > 0:
                self.lgfr = [self.rings[0]]
            else:
                self.lgfr = None

    def to_conformer(self):
        # TODO
        pass

    def to_rdmol(self):
        # I cannot use smilestomol as I need to know how atom ids are mapped
        nodes = sorted(list(self.graph.nodes(data=True)), key=lambda x: x[0])
        index_dict_no = {}  # d[atom idx in the new rdmol] = original graph node
        index_dict_on = {}

        atom_number_list = []
        for i in range(len(nodes)):
            graph_node = nodes[i][0]  # this could be different form i!
            symbol = nodes[i][1]['symbol']
            z = Element(symbol).Z
            atom_number_list.append(z)
            # if i == 0:
            #     mol = Chem.MolFromSmarts("[#" + str(z) + "]")
            #     rwmol = Chem.RWMol(mol)
            # else:
            #     rwmol.AddAtom(Chem.Atom(z))
            index_dict_no[i] = graph_node
            index_dict_on[graph_node] = i

        adj = nx.convert.to_dict_of_dicts(self.graph)

        new_ac = np.zeros((len(nodes), len(nodes))).astype(int)  # atomic connectivity
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if index_dict_no[j] in adj[index_dict_no[i]].keys():
                    new_ac[i, j] = 1
                    new_ac[j, i] = 1

        ap = ACParser(ac=new_ac, atomnumberlist=atom_number_list, charge=0)
        mol, smiles = ap.parse(charged_fragments=False, force_single=False, sanitize=True, expliciths=True)
        return mol, smiles, index_dict_no, index_dict_on

    def get_rings(self, k='nconnect', threshold=1):
        edges = [(u, v) for (u, v, d) in self.ring_graph.edges(data=True) if d[k] >= threshold]
        subgraph = self.ring_graph.edge_subgraph(edges).copy()
        return [[self.rings[iring] for iring in c] for c in
                sorted(nx.connected_components(subgraph), key=len, reverse=True)]

    @staticmethod
    def mol_to_nx(mol):
        """
        from https://github.com/maxhodak/keras-molecules/
        indices of nodes are inherited from atom indices in rdmol

        :param mol:
        :return:
        """
        G = nx.Graph()
        for atom in mol.GetAtoms():
            G.add_node(atom.GetIdx(),
                       # atomic_num=atom.GetAtomicNum(),
                       symbol=atom.GetSymbol(),
                       # formal_charge=atom.GetFormalCharge(),
                       # chiral_tag=atom.GetChiralTag(),
                       # hybridization=atom.GetHybridization(),
                       # num_explicit_hs=atom.GetNumExplicitHs(),
                       # is_aromatic=atom.GetIsAromatic(),
                       )
        for bond in mol.GetBonds():
            G.add_edge(
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                # bond_type=bond.GetBondType()
            )
        return G

    @classmethod
    def from_rdmol(cls, rdmol):
        smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)
        g = GraphMolecule.mol_to_nx(rdmol)
        if not nx.is_connected(g):
            warnings.warn('W: the graph is not connected!')
        return cls(g, smiles)

    def as_dict(self):
        # TODO
        d = OrderedDict()
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['graph_dict_of_lists'] = nx.to_dict_of_lists(self.graph)
        d['graph_node_symbols'] = OrderedDict(self.symbols)
        d['hsmiles'] = self.hsmiles
        return d

    def fingerprint(self):
        mol = Chem.MolFromSmiles(self.hsmiles)
        return FingerprintMols.FingerprintMol(mol)

    def fp_similarity(self, other, metric='Tanimoto'):
        """
        use RDK fingerprint similarity based on different metrics

        TODO add args to customize RDKfp, see https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint

        see Landrum2012 for more details

        :param str metric: "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet",
        :return:
        """
        mol = Chem.MolFromSmiles(self.hsmiles)
        other_mol = Chem.MolFromSmiles(other.hsmiles)
        fps = [FingerprintMols.FingerprintMol(mol), FingerprintMols.FingerprintMol(other_mol)]
        for func in DataStructs.similarityFunctions:
            if func[0] == metric:
                metric_function = func[1]
                return DataStructs.FingerprintSimilarity(fps[0], fps[1], metric=metric_function)
        return None

    @property
    def smarts(self):
        m = Chem.MolFromSmiles(self.hsmiles)
        return Chem.MolToSmarts(m)

    @property
    def n_nonsignlebond_electrons(self):
        m = Chem.MolFromSmiles(self.hsmiles)
        nve = Descriptors.NumValenceElectrons(m)
        nre = Descriptors.NumRadicalElectrons(m)
        nbonds = len(m.GetBonds())
        return nve - nre - 2 * nbonds

    @property
    def molweight(self):
        m = Chem.MolFromSmiles(self.hsmiles)
        return Descriptors.MolWt(m)

    def match(self, other):  # use smarts as patt
        if len(self.graph) == 1:
            nodes = list(self.graph.nodes(data=True))
            element = nodes[0][1]['symbol']
            if element in other.symbols.values():
                return True
        elif len(self.graph) == 2:
            nodes = list(self.graph.nodes(data=True))
            elements = {nodes[0][1]['symbol'], nodes[1][1]['symbol']}
            other_edges = list(other.graph.edges)
            other_elements = [{other.symbols[e[0]], other.symbols[e[1]]} for e in other_edges]
            if elements in other_elements:
                return True
        patt = Chem.MolFromSmarts(self.smarts)
        m = Chem.MolFromSmiles(other.hsmiles)
        return m.HasSubstructMatch(patt)

    @classmethod
    def from_dict(cls, d):
        g = nx.from_dict_of_lists(d['graph_dict_of_lists'])
        nx.set_node_attributes(g, d['graph_node_symbols'])
        return cls(g, d['hsmiles'])

    def __eq__(self, other):
        """
        #TODO test robustness

        :param other:
        :return:
        """
        g1 = self.graph
        g2 = other.graph
        return nx.is_isomorphic(g1, g2, node_match=iso.generic_node_match('symbol', None, eq))

    def __repr__(self):
        return self.hsmiles

    def partition(self, bone_selection='lgfr'):
        """
        parition the molecule into a backbone graph and a list of fragments (graphs)

        :param bone_selection:
        :return:
        """
        if bone_selection == 'lgfr':
            if self.lgfr is None:
                return None
            backbone_nodes = set([item for sublist in self.lgfr for item in sublist])
        else:
            if self.lgcr is None:
                return None
            backbone_nodes = set([item for sublist in self.lgcr for item in sublist])

        bone_graph = nx.Graph()
        fragments_graph = nx.Graph()  # this is the graph of anything other than bone nodes
        for node in self.graph:
            s = self.symbols[node]
            if node in backbone_nodes:
                bone_graph.add_node(node, symbol=s)

            else:
                fragments_graph.add_node(node, symbol=s)

        joints_bone_as_keys = {}
        joints_frag_as_keys = {}
        for n, nbrs in self.graph.adj.items():
            for nbr, d in nbrs.items():
                if n in backbone_nodes and nbr in backbone_nodes:
                    bone_graph.add_edge(n, nbr, **d)
                elif n not in backbone_nodes and nbr not in backbone_nodes:
                    fragments_graph.add_edge(n, nbr, **d)
                elif n in backbone_nodes and nbr not in backbone_nodes:
                    if n not in joints_bone_as_keys.keys():
                        joints_bone_as_keys[n] = [nbr]
                    else:
                        joints_bone_as_keys[n].append(nbr)
                elif n not in backbone_nodes and nbr in backbone_nodes:
                    if n not in joints_frag_as_keys.keys():
                        joints_frag_as_keys[n] = [nbr]
                    else:
                        joints_frag_as_keys[n].append(nbr)

        bone_graph.graph['joints'] = joints_bone_as_keys

        fragments = []
        for c in sorted(nx.connected_components(fragments_graph), key=len, reverse=True):
            nodes = list(c)
            joints_in_frag = {}
            for j in joints_frag_as_keys.keys():
                if j in nodes:
                    joints_in_frag[j] = joints_frag_as_keys[j]
            frag_graph = nx.Graph(joints=joints_in_frag)  # joints_in_frag[frag_node] is a list of bone joints
            frag_graph.add_nodes_from((n, self.graph.nodes[n]) for n in nodes)
            frag_graph.add_edges_from(
                (n, nbr, d) for n, nbrs in self.graph.adj.items() if n in nodes for nbr, d in nbrs.items() if
                nbr in nodes)
            fragments.append(frag_graph)
        return bone_graph, fragments

    def partition_to_bone_frags(self, bone_selection='lgfr'):
        bone_graph, fragments = self.partition(bone_selection)
        gb = GraphBackbone(bone_graph, self, partition_scheme=bone_selection)

        scs = [GraphSidechain(frag_graph, self) for frag_graph in fragments]
        scs = sorted(scs, key=lambda x: len(x.graph), reverse=True)
        return gb, scs


class GraphFragment:

    def __init__(self, subgraph, molgraph):
        self.graph = subgraph
        # self.molgraph = copy.deepcopy(molgraph)  # parent
        self.molgraph = molgraph  # parent
        self.joints = subgraph.graph['joints']
        self.rdmol, self.smarts = self.to_rdmol()
        self.symbols = nx.get_node_attributes(self.graph, 'symbol')

    def to_conformer(self):
        # TODO
        pass

    def __eq__(self, other):
        """
        #TODO test robustness

        :param other:
        :return:
        """
        g1 = self.graph
        g2 = other.graph
        return nx.is_isomorphic(g1, g2, node_match=iso.generic_node_match('symbol', None, eq))

    @property
    def n_nonsignlebond_electrons(self):
        nve = Descriptors.NumValenceElectrons(self.rdmol)
        nre = Descriptors.NumRadicalElectrons(self.rdmol)
        nbonds = len(self.rdmol.GetBonds())
        return nve - nre - 2 * nbonds

    @property
    def molweight(self):
        return Descriptors.MolWt(self.rdmol)

    def match_mol(self, molgraph):
        if len(self.graph) == 1:
            nodes = list(self.graph.nodes(data=True))
            element = nodes[0][1]['symbol']
            if element in molgraph.symbols.values():
                return True
        elif len(self.graph) == 2:
            nodes = list(self.graph.nodes(data=True))
            elements = {nodes[0][1]['symbol'], nodes[1][1]['symbol']}
            other_edges = list(molgraph.graph.edges)
            other_elements = [{molgraph.symbols[e[0]], molgraph.symbols[e[1]]} for e in other_edges]
            if elements in other_elements:
                return True
        patt = Chem.MolFromSmarts(self.smarts)
        m = Chem.MolFromSmiles(molgraph.hsmiles)
        return m.HasSubstructMatch(patt)

    def match_frag(self, frag):
        if len(self.graph) == 1:
            nodes = list(self.graph.nodes(data=True))
            element = nodes[0][1]['symbol']
            if element in frag.symbols.values():
                return True
        elif len(self.graph) == 2:
            nodes = list(self.graph.nodes(data=True))
            elements = {nodes[0][1]['symbol'], nodes[1][1]['symbol']}
            other_edges = list(frag.graph.edges)
            other_elements = [{frag.symbols[e[0]], frag.symbols[e[1]]} for e in other_edges]
            if elements in other_elements:
                return True
        patt = self.smarts
        m = Chem.MolFromSmarts(frag.smarts)
        return m.HasSubstructMatch(patt)

    def fp_similarity(self, other, metric='Tanimoto'):
        """
        use RDK fingerprint similarity based on different metrics

        TODO add args to customize RDKfp, see https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint

        see Landrum2012 for more details

        :param str metric: "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet",
        :return:
        """
        mol = self.rdmol
        other_mol = other.rdmol
        fps = [FingerprintMols.FingerprintMol(mol), FingerprintMols.FingerprintMol(other_mol)]
        for func in DataStructs.similarityFunctions:
            if func[0] == metric:
                metric_function = func[1]
                return DataStructs.FingerprintSimilarity(fps[0], fps[1], metric=metric_function)
        return None

    def to_rdmol(self, canonical=True):
        parent_rdmol, parent_smiles, parent_index_dict_no, parent_index_dict_on = self.molgraph.to_rdmol()

        nodes = sorted(list(self.graph.nodes(data=True)), key=lambda x: x[0])
        atom_number_list = []
        index_dict_no = {}  # d[atom idx in the new rdmol] = original graph node
        index_dict_on = {}

        for i in range(len(nodes)):
            graph_node = nodes[i][0]  # this is the atom idx in parent, could be different form i!
            symbol = nodes[i][1]['symbol']
            z = Element(symbol).Z
            atom_number_list.append(z)
            if i == 0:
                mol = Chem.MolFromSmarts("[#" + str(z) + "]")
                rwmol = Chem.RWMol(mol)
            else:
                rwmol.AddAtom(Chem.Atom(z))
            index_dict_no[i] = graph_node
            index_dict_on[graph_node] = i
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                atom_idx_i = index_dict_no[i]
                atom_idx_j = index_dict_no[j]
                if self.graph.has_edge(atom_idx_i, atom_idx_j):
                    bond = parent_rdmol.GetBondBetweenAtoms(atom_idx_i, atom_idx_j)
                    rwmol.AddBond(i, j, bond.GetBondType())

        mol = rwmol.GetMol()
        for i in range(len(nodes)):
            if nodes[i][0] in self.joints.keys():
                bone_joint_atom_idx = index_dict_no[i]
                side_joint_atom_ids = self.joints[bone_joint_atom_idx]
                a = mol.GetAtomWithIdx(i)
                radical_electrons_on_a = 0
                for side_joint_atom_idx in side_joint_atom_ids:
                    bond = parent_rdmol.GetBondBetweenAtoms(bone_joint_atom_idx, side_joint_atom_idx)
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        radical_electrons_on_a += 1
                    elif bond.GetBondType() == Chem.BondType.DOUBLE:
                        radical_electrons_on_a += 2
                    elif bond.GetBondType() == Chem.BondType.TRIPLE:
                        radical_electrons_on_a += 3
                    elif bond.GetBondType() == Chem.BondType.AROMATIC:
                        radical_electrons_on_a += 2
                a.SetNumRadicalElectrons(radical_electrons_on_a)

        smarts = Chem.MolToSmarts(mol)
        if canonical:
            expliciths = True
            sanitize = True
            smiles = Chem.MolToSmiles(mol, allHsExplicit=expliciths, isomericSmiles=True)
            mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
            smarts = Chem.MolToSmarts(mol)
        return mol, smarts

    def fingerprint(self):
        return FingerprintMols.FingerprintMol(self.rdmol)

    def as_dict(self):
        # TODO
        d = OrderedDict()
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['graph_dict_of_lists'] = nx.to_dict_of_lists(self.graph)
        d['graph_node_symbols'] = OrderedDict(self.symbols)
        d['parent_graph'] = self.molgraph.to_dict()
        d['smarts'] = self.smarts
        d['joints'] = self.joints
        return d

    @classmethod
    def from_dict(cls, d):
        g = nx.from_dict_of_lists(d['graph_dict_of_lists'])
        nx.set_node_attributes(g, d['graph_node_symbols'])
        return cls(g, GraphMolecule.from_dict(d['parent_graph']))


class GraphBackbone(GraphFragment):
    def __init__(self, subgraph, molgraph, partition_scheme='lgcr'):
        super().__init__(subgraph, molgraph)
        self.unsaturation = self.n_nonsignlebond_electrons / len(self.graph)
        # if partition_scheme == 'lgcr':
        #     self.rings = molgraph.lgcr
        # elif partition_scheme == 'lgfr':
        #     self.rings = molgraph.lgfr
        self.rings = nx.minimum_cycle_basis(self.graph)  # this should be identical to lgcr or lgfr
        self.rings = sorted(self.rings, key=lambda x: len(x))
        self.nrings = len(self.rings)
        self.partition_scheme = partition_scheme

    def as_dict(self):
        # TODO
        d = OrderedDict()
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['graph_dict_of_lists'] = nx.to_dict_of_lists(self.graph)
        d['graph_node_symbols'] = OrderedDict(self.symbols)
        d['parent_graph'] = self.molgraph.to_dict()
        d['smarts'] = self.smarts
        d['joints'] = self.joints
        d['nrings'] = self.nrings
        d['partition_scheme'] = self.partition_scheme
        d['unsaturation'] = self.unsaturation
        d['rings'] = self.rings
        return d

    @classmethod
    def from_dict(cls, d):
        g = nx.from_dict_of_lists(d['graph_dict_of_lists'])
        nx.set_node_attributes(g, d['graph_node_symbols'])
        return cls(g, GraphMolecule.from_dict(d['parent_graph']), d['parition_scheme'])


class GraphSidechain(GraphFragment):
    def __init__(self, subgraph, molgraph):
        """
        it is not possible for a side chain to have two frag-joints connected to one bone-joint, as this should be
        already in backbone

        side chain (frags) cannot have 2 side joints, each connected to one bone-joint as this would be a ring fused
         to backbone (side chain is a connected component)
        :var self.joints: a dict, keys are nodes in the subgraph that have parent graph edges not in this subgraph

        """
        super().__init__(subgraph, molgraph)
        self.unsaturation = self.n_nonsignlebond_electrons / len(self.graph)
        if len(self.joints.keys()) > 1:
            warnings.warn('W: sidechain has more than one side joint... impossible!')
        sc_joint = list(self.joints.keys())[0]
        for node in subgraph:
            rank = nx.shortest_path_length(subgraph, source=node, target=sc_joint)
            subgraph.nodes[node]['rank'] = rank
        self.rings = nx.minimum_cycle_basis(self.graph)  # technically sssr
        self.rings = sorted(self.rings, key=lambda x: len(x))
        self.nrings = len(self.rings)
        if len(self.rings) > 0:
            self.hasring = True
        else:
            self.hasring = False

    def as_dict(self):
        # TODO
        d = OrderedDict()
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        d['graph_dict_of_lists'] = nx.to_dict_of_lists(self.graph)
        d['graph_node_symbols'] = OrderedDict(self.symbols)
        d['parent_graph'] = self.molgraph.to_dict()
        d['smarts'] = self.smarts
        d['joints'] = self.joints
        d['nrings'] = self.nrings
        d['unsaturation'] = self.unsaturation
        d['rings'] = self.rings
        return d

    @classmethod
    def from_dict(cls, d):
        g = nx.from_dict_of_lists(d['graph_dict_of_lists'])
        nx.set_node_attributes(g, d['graph_node_symbols'])
        return cls(g, GraphMolecule.from_dict(d['parent_graph']))
