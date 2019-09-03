from pymatgen.core.structure import *
from ELEMENTS import ELEMENTS
import numpy as np
import pickle
import sys
import anytree
import copy

def get_pbc(infilename):
    lattice = Structure.from_file(infilename).lattice
    a = lattice.matrix[0]
    b = lattice.matrix[1]
    c = lattice.matrix[2]
    return np.array([a,b,c])

def get_nbs(site, mol, co=1.5):
    """
    get a list of nbs of a site in a certain mol
    :param site:
    :param mol:
    :param co: coefficient of determining cutoff
    :return: a list of nbs of a site in a certain mol
    """
    all_sites = mol.sites
    nbs= []
    for i in all_sites:
        if not np.allclose(i.coords, site.coords):
            cutoff = (ELEMENTS[str(site.specie)][1] + ELEMENTS[str(i.specie)][1]) * co
            distance = np.linalg.norm(i.coords - site.coords)
            if distance < cutoff:
                nbs.append(i)
    return nbs

def write_obj(obj, objnamestring):
    with open(objnamestring, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(obj, f)

def read_obj(pklfilename):
    with open(pklfilename, 'rb') as f:  # Python 3: open(..., 'rb')
        obj = pickle.load(f)
    return obj

def unify(vector):
    # print (inspect.stack()[1][3])
    if np.linalg.norm(vector) == 0.0:
        print('v:', vector)
        sys.exit(1)
    else:
        return vector/np.linalg.norm(vector)

def angle_btw(v1, v2):
    v1_u = unify(v1)
    v2_u = unify(v2)
    angle_radian = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return angle_radian

def analyze_bone(filename):
    # s1
    # |
    # s2
    mol = Molecule.from_file(filename)

    ends = []
    for site in mol.sites:
        nbs = get_nbs(site, mol)
        if len(nbs) == 2 and site.species_string == 'C':
            ends.append(site)
    s1, s2 = ends
    outv_s1 = unify(s1.coords - s2.coords)
    outv_s2 = -1 * outv_s1

    labeled_sites = []
    for site in mol.sites:
        if np.allclose(site.coords, s1.coords):
            ns = Site(site.species_string, site.coords, properties={'label':'C_cjunction1'})
        elif np.allclose(site.coords, s2.coords):
            ns = Site(site.species_string, site.coords, properties={'label': 'C_cjunction2'})
        else:
            ns = Site(site.species_string, site.coords, properties={'label': 'core'})
        labeled_sites.append(ns)
    labeled_mol = Molecule.from_sites(labeled_sites)

    return [s1, s2, outv_s1, outv_s2, labeled_mol]

def vtree2bltree(cj_node):
    for level in reversed(list(anytree.LevelOrderGroupIter(cj_node))):
        for node in level:
            if not node.is_root:
                node.coords = np.linalg.norm(node.coords)
    return cj_node

def rotation_matrix(axis, theta_d):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    np.dot(rotation_matrix(axis,theta_d), v)
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    import math
    theta = np.deg2rad(theta_d)
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                    [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                    [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotate_along_axis(v, axis, theta_d):
    return np.dot(rotation_matrix(axis,theta_d), v)

def append_sidechain_number(node, number):
    nnode = copy.deepcopy(node)
    for n in list(anytree.LevelOrderIter(nnode))[1:]:
        n.symbol = n.symbol + '_' + number
    return nnode

def align_decompose_bltree2vtree(cj_node, v_extending):
    #
    # p --> n --> c1, c2, c3
    #

    cj_n = copy.deepcopy(cj_node)
    cj_n.children[0].coords = cj_n.children[0].coords * v_extending
    sj_node = cj_n.children[0]
    for n in list(anytree.LevelOrderIter(sj_node)):
        v_p2n = unify(n.coords)
        if len(n.children) == 1:
            n.children[0].coords =  n.children[0].coords * v_p2n
        elif len(n.children) == 3:
            trivial_v1 = [1,0,0]
            trivial_v2 = [0,1,0]
            if angle_btw(v_p2n, trivial_v1) != 0 and angle_btw(v_p2n, trivial_v1) != 180:
                v_aux = unify(np.cross(trivial_v1, v_p2n))
            else:
                v_aux = unify(np.cross(trivial_v2, v_p2n))
            v_n2c1 = rotate_along_axis(v_p2n, v_aux, 180.0-109.5)
            v_n2c2 = rotate_along_axis(v_n2c1, v_p2n, 120.0)
            v_n2c3 = rotate_along_axis(v_n2c2, v_p2n, 120.0)
            n.children[0].coords = n.children[0].coords * v_n2c1
            n.children[1].coords = n.children[1].coords * v_n2c2
            n.children[2].coords = n.children[2].coords * v_n2c3
        elif len(n.children) != 0:
            print('ERROR: in the side tree there are nodes other than sp3 or sp')
            sys.out(1)
    return cj_n

def rotateset_3vs(vs, ang_res): # rotate 3 vectors
    axis = vs[0] + vs[1] + vs[2]
    theta_d_set = np.arange(0, 360, ang_res)
    new_vs_set = []
    for theta_d in theta_d_set:
        vn0 = rotate_along_axis(vs[0], axis, theta_d)
        vn1 = rotate_along_axis(vs[1], axis, theta_d)
        vn2 = rotate_along_axis(vs[2], axis, theta_d)
        new_vs_set.append([vn0, vn1, vn2])
    return new_vs_set

def adds2rot3(in_str, ang_res, node_to_be_populated):
    in_sites = in_str.sites
    vs = [node_to_be_populated.children[0].coords, node_to_be_populated.children[1].coords, node_to_be_populated.children[2].coords]
    for site in in_sites:
        if 'label' in site.properties.keys():
            if site.properties['label'] == node_to_be_populated.symbol:
                current_site = site
    new_vs_set = rotateset_3vs(vs, ang_res)
    set_new_str = []
    ref_coords = current_site.coords
    ns1_str = node_to_be_populated.children[0].symbol.split('_')[0]
    ns2_str = node_to_be_populated.children[1].symbol.split('_')[0]
    ns3_str = node_to_be_populated.children[2].symbol.split('_')[0]
    for v in new_vs_set:
        ns1 = Site(ns1_str, v[0] + ref_coords, properties={'label':node_to_be_populated.children[0].symbol})
        ns2 = Site(ns2_str, v[1] + ref_coords, properties={'label':node_to_be_populated.children[1].symbol})
        ns3 = Site(ns3_str, v[2] + ref_coords, properties={'label':node_to_be_populated.children[2].symbol})
        adding_sites = [ns1, ns2, ns3]
        new_sites = in_sites + adding_sites
        set_new_str.append(Molecule.from_sites(new_sites))
    return set_new_str

def adds2rigid(in_str, ang_res, node_to_be_populated):
    in_sites = in_str.sites
    for site in in_sites:
        if 'label' in site.properties.keys():
            if site.properties['label'] == node_to_be_populated.symbol:
                current_site = site
    #print(current_site, current_site.properties)
    ref_coords = current_site.coords
    adding_sites = []
    for child in node_to_be_populated.children:
        ns = Site(child.symbol.split('_')[0], child.coords + ref_coords, properties={'label':child.symbol})
        adding_sites.append(ns)
    set_new_str = [Molecule.from_sites(adding_sites + in_sites)]
    return set_new_str

def translate_a_copy(mol,transv):
    new_copy = []
    for site in mol.sites:
        new_site = Site(site.species_string, site.coords + transv)
        new_copy.append(new_site)
    return Molecule.from_sites(new_copy)

def copies_too_close(mol, trans_vs):
    too_close_yon = False
    for v in trans_vs:
        copy = translate_a_copy(mol, v)
        if too_close(mol, copy):
            too_close_yon = True
            break
    return too_close_yon


def survive(mollist, pbc):
    '''
    :param mollist:
    :param pbc: a list of vectors in cart
    :return:
    '''
    a, b, c = pbc
    supercell_trans_v = [a, b, c, a+b, a+c, a+b+c, b+c]
    survival = []
    for mol in mollist:
        if not self_too_close(mol):
            if not copies_too_close(mol, supercell_trans_v):
                survival.append(mol)
    return survival

def get_discutoff(s1, s2):
    dis_cutoff = (ELEMENTS[str(s1.specie)][1] + ELEMENTS[str(s2.specie)][1]) * 1.5
    return dis_cutoff

def self_too_close(mol):
    too_close = False
    for site in mol.sites:
        nbs = get_nbs(site, mol)
        if site.species_string == 'C' and len(nbs) > 4:
            too_close = True
            break
        elif site.species_string == 'Ge' and len(nbs) > 4:
            too_close = True
            break
        elif site.species_string == 'Si' and len(nbs) > 4:
            too_close = True
            break
        elif site.species_string == 'H' and len(nbs) > 2:
            too_close = True
            break
    return too_close


def too_close(mol1, mol2):
    yon = False
    for s1 in mol1.sites:
        for s2 in mol2.sites:
            if s1.distance(s2) < get_discutoff(s1, s2):
                yon = True
                return yon
    return yon

def grow_test(in_str, core_junction_node):
    grow_t=[]
    for node in anytree.LevelOrderIter(core_junction_node):
        print('we are at')
        print(node)
        if len(node.children) < 3:
            set_new_str = adds2rigid(in_str, 0, node)
        else:
            set_new_str = adds2rot3(in_str, 360, node)
        in_str = set_new_str[0]
        grow_t.append(set_new_str[0])
    for i in range(len(grow_t)):
        grow_t[i].to('xyz', 'grow_' + str(i) + '.xyz')
    return


from timeit import default_timer as timer


def grow(in_set, ang_res, cj_node1, cj_node2, pbc, mol_sym):
    if mol_sym == 'inversion':
        node1_list = list(anytree.LevelOrderGroupIter(cj_node1))
        node2_list = list(anytree.LevelOrderGroupIter(cj_node2))
        node_list_combine = []
        for i in range(len(node1_list)):
            group = node1_list[i] + node2_list[i]
            node_list_combine.append(group)
    for level in node_list_combine:
        for node in level:
            start = timer()
            survival = []
            for structure in in_set:
                if len(node.children) == 3:
                    if node.children[0].symbol.split('_')[0] == 'H' and \
                        node.children[1].symbol.split('_')[0] == 'H' and \
                        node.children[2].symbol.split('_')[0] == 'H':
                            set_new_str = adds2rigid(structure, ang_res, node)
                    else:
                        set_new_str = adds2rot3(structure, ang_res, node)
                elif len(node.children) <3:
                    set_new_str = adds2rigid(structure, ang_res, node)
                #print(set_new_str[0])
                survival = survival + survive(set_new_str, pbc)
                #print(survival[0])
            in_set = survival
            print('we are at')
            print(node)
            print('survival:',len(survival))
            end = timer()
            print(end - start)
        #print(anytree.RenderTree(node))
    return survival

def pack_mol2pbc(mol, pbc):
    species = []
    coords = []
    pbcsites = []
    for site in mol.sites:
        species.append(site.species_string)
        coords.append(site.coords)
    str = Structure(pbc, species, coords, coords_are_cartesian=True)
    return str