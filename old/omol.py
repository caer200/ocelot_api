import copy
import itertools
import anytree
import mathop
import sitesop
from sidechain import SideChain
from backbone import Backbone
from msite import MSite
from scipy.spatial.distance import pdist, squareform

class OMol:
    """
    OMol  --  sidechain1
          --  sidechain2
          --  ...
          --  backbone
    """

    def __init__(self, sites):
        """
        after init each site should have...
        site.mid                                        an int starts from 0
        site.properties['num_neighbors']                an int
        site.properties['hybrid']                       'sp3', 'sp2', 'sp', 'ion', 'isolated', 'hydrogen'
        site.properties['relative_position']            'bone', 'side', 'bone-joint', 'side-joint'
        site.properties['sidechain_label']              [None, None], [int, int] both start from 0,
                                                        [sidechain_counter, level], level is 0 for bone-joint

        we will use 'relative_position' and 'sidechain_label' to attach generic sidechains to backbone in polygen
        the vector between bone-joint and side-joint will determine all other bond angles
        :param sites:
        """

        self.sites = tuple(sites)
        self.coordmat = sitesop.coordmat(sites)
        self.distmat = squareform(pdist(self.coordmat))
        self.bondmat = sitesop.bondmat(self.sites, self.distmat)
        self.nrnbrmap = sitesop.nrnbrmap(self.sites)
        self.nsites = len(self.sites)

        for i in range(self.nsites):
            self.sites[i].mid = i
            self.sites[i].properties['num_neighbors'] = len(self.nrnbrmap[i])
            self.sites[i].properties['hybrid'] = OMol.check_hybid(self.sites[i].symbol, len(self.nrnbrmap[i]))

        self.side_idx = []
        self.bone_idx = []
        # partition mol into backbone and side chains based on rings
        # TODO the group ring method cannot handle rubrene
        rings = mathop.loop_finding(self.nrnbrmap, 6) + mathop.loop_finding(self.nrnbrmap, 5) + \
                mathop.loop_finding(self.nrnbrmap, 4)
        ring_site_idx = [i for ring in rings for i in ring]
        largest_group = sorted(self.group_bonded_sites(self.bondmat, ring_site_idx), key=lambda x: len(x))[-1]

        for i in range(len(self.sites)):
            # TODO not sound here, need a better way to identify backbone
            if len([idx for idx in self.nrnbrmap[i] if
                    self.sites[idx].properties['hybrid'] in ['sp2'] or self.sites[idx].symbol not in ['C', 'H']]) > 1\
                    and i in largest_group:
                self.sites[i].properties['relative_position'] = 'bone'
                self.sites[i].properties['sidechain_label'] = [None, None]
                self.bone_idx.append(i)
            else:
                self.sites[i].properties['relative_position'] = 'side'
                self.sites[i].properties['sidechain_label'] = [-1, -1]
                self.side_idx.append(i)

        # here we assume there's only one atom extended out from each backbone joint
        # TODO 68.cif  this is not the case
        # the property sidechain_label denotes [the label of the side chain, the level of that site on side chain]
        sidechain_counter = 0
        for i in range(len(self.sites)):
            if not set(self.nrnbrmap[i]).issubset(set(self.bone_idx)) and i in self.bone_idx:
                self.sites[i].properties['relative_position'] = 'bone-joint'
                for j in self.nrnbrmap[i]:
                    if j not in self.bone_idx:
                        self.sites[j].properties['relative_position'] = 'side-joint'
                        self.sites[i].properties['sidechain_label'] = [sidechain_counter, 0]
                        self.sites[j].properties['sidechain_label'] = [sidechain_counter, 1]
                        # break
                sidechain_counter += 1
        self.num_sidechains = sidechain_counter

        # sitesop.sites_toxyz([self.msites[i] for i in self.bone_idx], 'bone.xyz')

        # now assign sidechain_label to non-joint msites on side chains
        bonded_side_sites = OMol.group_bonded_sites(self.bondmat, self.side_idx)
        idxtrees = [None] * self.num_sidechains
        for grp in bonded_side_sites:
            grp_sidechain_counter = -1

            # print(grp)
            # for i in grp:
            #     print(self.msites[i])

            for idx in grp:
                if self.sites[idx].properties['sidechain_label'][0] != -1 and \
                        self.sites[idx].properties['sidechain_label'][1] == 1:
                    grp_sidechain_counter = self.sites[idx].properties['sidechain_label'][0]
                    break

            bone_joint = [i for i in self.bone_idx if
                          self.sites[i].properties['relative_position'] == 'bone-joint'
                          and self.sites[i].properties['sidechain_label'][0] == grp_sidechain_counter][0]

            for idx in grp:
                self.sites[idx].properties['sidechain_label'][0] = grp_sidechain_counter

            # now assign the 2nd inx of sidechain label
            bone_joint_node = OMol.chain2sidetree(bone_joint, self.bondmat, grp)
            # this is redundant as levels have been calculated in chain2sidetree
            levels = [[node.name for node in children] for children in anytree.LevelOrderGroupIter(bone_joint_node)]
            for ilevel in range(len(levels)):
                for idx in levels[ilevel]:
                    self.sites[idx].properties['sidechain_label'][1] = ilevel
            idxtrees[grp_sidechain_counter] = bone_joint_node


        self.sitetrees = [OMol.idxtree2sitetree(tree, self) for tree in idxtrees]

        self.backbone = Backbone(self, [self.sites[i] for i in self.bone_idx])

        self.sidechains = [SideChain(self, tree) for tree in self.sitetrees]  # this muse be after backbone!

        self.data = {
            'backbone' : self.backbone.data,
            'sidechains' : [sc.data for sc in self.sidechains],
        }

    def to_xyz(self, xyzname):
        sitesop.sites_toxyz(self.sites, xyzname)

    def __repr__(self):
        outs = ['OMol']
        for s in self.sites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @classmethod
    def from_xyz(cls, xyz):
        with open(xyz, 'r') as f:
            lines = f.readlines()[2:]
        sites = []
        for line in lines:
            items = line.strip().split()
            if len(items) == 4:
                s = MSite(items[0], items[1:])
                sites.append(s)
        return cls(sites)

    @staticmethod
    def idxtree2sitetree(idxtree, mol):
        symboltree = copy.deepcopy(idxtree)
        for node in anytree.LevelOrderIter(symboltree):
            node.site = mol.sites[node.name]
        return symboltree

    @staticmethod
    def chain2sidetree(root_idx, bondmat, idxrange):
        """
        root_idx should only have one bonded idx in idxrange, this is the grow direction
        notice nodes are mutable
        # TODO add compatibility of closed loop
        :param root_idx: bone-joint
        :param bondmat:
        :param idxrange: the idx for one group, this does not include root_idx which is on the bone
        :return:
        """
        idxrange = [root_idx] + idxrange  # root_idx must come first
        # initialization
        classified = [[root_idx], ]
        pointer = 0
        while pointer != len(classified):
            outside = [idx for idx in idxrange if idx not in list(itertools.chain.from_iterable(classified))]
            tobeincluded = []
            for j in classified[pointer]:
                for i in outside:
                    if bondmat[j][i]:
                        tobeincluded.append(i)
            if len(tobeincluded) != 0:
                classified += [tobeincluded]
            pointer += 1

        classified_nodes = []
        for i in classified:
            level = []
            for j in i:
                level.append(anytree.Node(j))
            classified_nodes.append(level)

        # now convert to tree
        for i in reversed(range(1, len(classified_nodes))):
            for j in range(len(classified_nodes[i])):
                upper = i - 1
                for k in range(len(classified_nodes[upper])):
                    if bondmat[classified_nodes[i][j].name][classified_nodes[upper][k].name]:
                        classified_nodes[i][j].parent = classified_nodes[upper][k]
        return classified_nodes[0][0]

    @staticmethod
    def group_bonded_sites(bondmat, idxrange):
        """
        :param bondmat: bondmat[i][j] returns whether ij are bonded
        :param idxrange:
        :return:
        """
        visited = []
        block_list = []
        while len(visited) != len(idxrange):
            # initialization
            unvisited = [idx for idx in idxrange if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in idxrange if idx not in block and idx not in visited]
                for i in outside:
                    if bondmat[block[pointer]][i]:
                        block.append(i)
                visited.append(block[pointer])
                pointer += 1
            block_list.append(block)
        return block_list

    @staticmethod
    def isomol():
        # TODO identify whether a molecule can be considered as an OMol
        pass

    @staticmethod
    def check_hybid(symbol, nnbs):
        if symbol == 'H':
            return 'hydrogen'
        elif symbol in ['O', 'S', 'Se']:
            if nnbs == 2:
                return 'sp3'
            elif nnbs == 1:
                return 'sp2'
            elif nnbs == 0:
                return 'isolated'
        elif symbol in ['N', 'P', 'As']:
            if nnbs == 3:
                return 'sp3'
            elif nnbs == 2:
                return 'sp2'
            elif nnbs == 1:
                return 'sp'
            elif nnbs == 0:
                return 'isolated'
        elif symbol in ['C', 'Si', 'Ge']:
            if nnbs == 4:
                return 'sp3'
            elif nnbs == 3:
                return 'sp2'
            elif nnbs == 2:
                return 'sp'
            elif nnbs == 1:
                return 'ion'
            elif nnbs == 0:
                return 'isolated'
