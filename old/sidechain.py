import anytree
import sitesop
import copy
import mathop
import numpy as np


class SideChain:
    """
    a generic scheme to generate the side chain
    # TODO add molecular surface calculation
    """

    def __init__(self, omol, sitetree):
        """
        :param omol: parent obj
        :param sitetree: bone-joint node, imported from omol
        """
        self.omol = omol
        self.raw_sitetree = sitetree
        self.sites = [node.site for node in anytree.PreOrderIter(self.raw_sitetree) if
                      node.site.properties['relative_position'] != 'bone-joint']
        self.volume = sitesop.volume(self.sites)
        self.sidechain_number = self.sites[0].properties['sidechain_label'][0]
        self.ishydrogen = len(self.sites) == 1 and self.sites[0].symbol == 'H'
        self.natoms = len(self.sites)
        self.nelectrons = sitesop.nelectrons(self.sites)

        # decide the angle between v_sidejoint-backbongom and v_backbone_p
        self.attached_bonejoint = next(s for s in self.omol.sites if s.properties['relative_position']
                                       == 'bone-joint' and s.properties['sidechain_label'][0] == self.sidechain_number)
        self.angle_vp = mathop.angle_btw(self.attached_bonejoint.coords - self.omol.backbone.gom,
                                         self.omol.backbone.vp, output='degree')
        if self.angle_vp > 90.0:
            self.angle_vp = 180 - self.angle_vp

        # let's look at morphology
        if self.natoms != 1:
            for node in anytree.LevelOrderIter(self.raw_sitetree):
                if len(node.children) > 1 or len(node.children) == 0:
                    self.branchnode = copy.deepcopy(node)  # like the silicon in tips
                    break
        else:
            self.branchnode = self.raw_sitetree.children[0]
        self.branchsites = [node.site for node in anytree.LevelOrderIter(self.branchnode)]  # add branch node
        self.branchvolume = sitesop.volume(self.branchsites)
        self.edensity = self.nelectrons / self.volume
        self.edensity_branch = self.nelectrons / self.branchvolume

        norms = np.zeros(len(self.branchsites))
        bs = self.branchnode.site
        for i in range(len(self.branchsites)):
            norms[i] = bs.distance(self.branchsites[i])
        self.branchnorms = list(norms)

        r = max(self.branchnorms)
        self.branchspherev = np.pi * r ** 3 * 4 / 3

        self.data = {
            'edensity' : self.edensity,
            'edensity_branch' : self.edensity_branch,
            'branchnorms' : self.branchnorms,
            'branchspherev' : self.branchspherev,
            'branchvolume' : self.branchvolume,
            'volume' : self.volume,
            'angle_vp' : self.angle_vp,
            'sidechain_number' : self.sidechain_number,
            'ishydrogen' : self.ishydrogen,
            'natoms' : self.natoms,
            'nelectrons' : self.nelectrons,
        }

