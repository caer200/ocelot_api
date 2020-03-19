import itertools
from collections import OrderedDict
from operator import eq
import pynauty as pn
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from networkx import Graph
from pymatgen.core.periodic_table import Element
from pymatgen.vis.structure_vtk import EL_COLORS
import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D as D3d
from rdkit.Chem.Fingerprints import FingerprintMols

"""
commonly used rdkit functions
"""

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
from rdkit.Chem import Mol
from rdkit.Chem import ResonanceMolSupplier
def makeRd(smile):
    rdmol = Chem.MolFromSmiles(smile)
    rdmol = Chem.AddHs(rdmol)
    return rdmol;

def HASH(graph):
        """
        see https://stackoverflow.com/questions/46999771/
        use with caution...

        :return:
        """
        g: nx.Graph = graph
        
        pnGraph=pn.Graph(g.number_of_nodes());
        edg=list(g.edges);
        nodesColored=[];
        colors={"H":[],"He":[],"Li":[],"Be":[],"B":[],"C":[],"N":[],"O":[],"F":[],"Ne":[],"Na":[],"Mg":[],"Al":[],"Si":[],"P":[],"S":[],"Cl":[],"Ar":[],"K":[],"Ca":[],"Sc":[],"Ti":[],"V":[],"Cr":[],"Mn":[],"Fe":[],"Co":[],"Ni":[],"Cu":[],"Zn":[],"Ga":[],"Ge":[],"As":[],"Se":[],"Br":[],"Kr":[],"Rb":[],"Sr":[],"Y":[],"Zr":[],"Nb":[],"Mo":[],"Tc":[],"Ru":[],"Rh":[],"Pd":[],"Ag":[],"Cd":[],"In":[],"Sn":[],"Sb":[],"Te":[],"I":[],"Xe":[],"Cs":[],"Ba":[],"La":[],"Ce":[],"Pr":[],"Nd":[],"Pm":[],"Sm":[],"Eu":[],"Gd":[],"Tb":[],"Dy":[],"Ho":[],"Er":[],"Tm":[],"Yb":[],"Lu":[],"Hf":[],"Ta":[],"W":[],"Re":[],"Os":[],"Ir":[],"Pt":[],"Au":[],"Hg":[],"Tl":[],"Pb":[],"Bi":[],"Po":[],"At":[],"Rn":[],"Fr":[],"Ra":[],"Ac":[],"Th":[],"Pa":[],"U":[],"Np":[],"Pu":[],"Am":[],"Cm":[],"Bk":[],"Cf":[],"Es":[],"Fm":[],"Md":[],"No":[],"Lr":[],"Rf":[],"Db":[],"Sg":[],"Bh":[],"Hs":[],"Mt":[],"Ds":[],"Rg":[],"Cn":[],"Nh":[],"Fl":[],"Mc":[],"Lv":[],"Ts":[],"Og":[]};
        for E in edg:
            pnGraph.connect_vertex(E[0],E[1]);
            try:
                nodesColored.index(E[0]);
                try:
                    colors[g.nodes[E[0]]["symbol"]].append(E[0]);
                except KeyError:
                    colors[g.nodes[E[0]]["symbol"]]=[];
                    colors[g.nodes[E[0]]["symbol"]].append(E[0]);
            except ValueError:
                nodesColored.append(E[0]);
                try:
                    colors[g.nodes[E[0]]["symbol"]].append(E[0]);
                except KeyError:
                    colors[g.nodes[E[0]]["symbol"]]=[];
                    colors[g.nodes[E[0]]["symbol"]].append(E[0]);
            try:
                nodesColored.index(E[1]);
                try:
                    colors[g.nodes[E[1]]["symbol"]].append(E[1]);
                except KeyError:
                    colors[g.nodes[E[1]]["symbol"]]=[];
                    colors[g.nodes[E[1]]["symbol"]].append(E[1]);
            except ValueError:
                nodesColored.append(E[1]);
                try:
                    colors[g.nodes[E[1]]["symbol"]].append(E[1]);
                except KeyError:
                    colors[g.nodes[E[1]]["symbol"]]=[];
                    colors[g.nodes[E[1]]["symbol"]].append(E[1]);
        j=-1;
        for c in colors:
            j=j+1;
            print(str(c)+" "+str(colors[c]))
            if colors[c]!=[]:
                pnGraph.set_vertex_coloring([set(colors[c])]);
            else:
                pnGraph.set_vertex_coloring([set([])]);
        return hash(pn.certificate(pnGraph));
    
def from_rdmol(rdmol, atomidx2nodename=None):
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
        raise Exception('the graph is not connected!')
    return g
print(HASH(from_rdmol(makeRd("CCF"))))
print(HASH(from_rdmol(makeRd("CCBr"))))