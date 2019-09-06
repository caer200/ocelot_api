"""
identify packing pattern
the input should be a pbc config of H terminated backbones

1. heuristic rules 10.1039/c7tc02553j, that worked only for small acenes, we want more

2. finger print Hirshfield 10.1039/c1ce05763d
"""

from collections import OrderedDict
# from schema.dimercollection import DimerCollection

class PackingIdentifier:

    def __init__(self, boneconfig):
        self.boneconfig = boneconfig

    def identify_heuristic(self):
        """
        :return:
        """
        bc_dimers, bc_dimers_transv_fcs = self.boneconfig.dimers_array, self.boneconfig.transv_fcs
        report = OrderedDict()
        for i in range(self.boneconfig.z):
            n_close_azm_and_parallel = 0  # min distance between backbone proj<2, 1.5>oslip
            n_close_azm_and_notparallel = 0  # min distance between backbone proj<2, 1.5>oslip
            n_close_vertical_and_parallel = 0  # min distance between backbone proj<2, 4>oslip>1.5
            n_close_vertical_and_notparallel = 0
            n_parallel_and_overlap = 0
            n_notparallel_and_overlap = 0

            dimers_ref_i = bc_dimers[i].flatten()

            dimers_ref_i = sorted(dimers_ref_i, key=lambda x: x.vslipnorm)[:27 * self.boneconfig.z]
            # DimerCollection(dimers_ref_i).to_xyz('dimers_{}.xyz'.format(i))

            for d in dimers_ref_i:
                if d.is_not_identical and d.is_close:
                    # if d.is_not_identical and d.is_close:
                    overlap, refboneproj, varboneproj, mindist2d = d.bone_overlap()
                    # print(d.minbonedist)
                    # print(mindist2d)
                    # debug += 1
                    # d.to_xyz('dimer_{}.xyz'.format(debug))
                    # overlap, refboneproj, varboneproj = d.mol_overlap()
                    if 1e-5 < mindist2d < 2.5:  # exclude overlap as dist=0.0 that case
                        if 4 > d.oslip > 1.5:
                            if d.oangle < 10.0:
                                n_close_vertical_and_parallel += 1
                            else:
                                n_close_vertical_and_notparallel += 1
                        elif d.oslip <= 1.5:
                            if d.oangle < 10.0:
                                n_close_azm_and_parallel += 1
                            else:
                                n_close_azm_and_notparallel += 1
                    if overlap > 1e-5:
                        if d.oangle < 10.0:
                            n_parallel_and_overlap += 1
                        else:
                            n_notparallel_and_overlap += 1

            if n_parallel_and_overlap > 2:
                packing = 'brickwork'

            elif n_close_vertical_and_parallel > 2 or n_close_azm_and_parallel > 2 or n_close_vertical_and_parallel + n_close_azm_and_parallel + n_parallel_and_overlap > 2:
                packing = 'edge_brickwork'

            elif n_parallel_and_overlap == 2:
                packing = 'slipped_stack'

            elif n_close_vertical_and_parallel == 2 or n_parallel_and_overlap + n_close_vertical_and_parallel + n_close_azm_and_parallel == 2:
                packing = 'edge_slipped_stack'

            elif n_notparallel_and_overlap >= 1:
                packing = 'herringbone'

            elif n_close_vertical_and_notparallel >= 1 or n_close_azm_and_notparallel >= 1:
                packing = 'edge_herringbone'

            elif n_parallel_and_overlap == 1 or n_close_vertical_and_parallel == 1:
                packing = 'dimerized'

            elif n_notparallel_and_overlap == n_parallel_and_overlap == n_close_azm_and_notparallel == n_close_azm_and_parallel == n_close_vertical_and_parallel == n_close_vertical_and_notparallel == 0:
                packing = 'sparse'

            else:
                packing = 'other'

            report[i] = OrderedDict(
                n_close_azm_and_parallel=n_close_azm_and_parallel,  # min distance between backbone proj<2, 1.5>oslip
                n_close_azm_and_notparallel=n_close_azm_and_notparallel,
                # min distance between backbone proj<2, 1.5>oslip
                n_close_vertical_and_parallel=n_close_vertical_and_parallel,
                # min distance between backbone proj<2, 4>oslip>1.5
                n_close_vertical_and_notparallel=n_close_vertical_and_notparallel,
                n_parallel_and_overlap=n_parallel_and_overlap,
                n_notparallel_and_overlap=n_notparallel_and_overlap,
                packing=packing,
            )

        return report
        # return [report[i]['packing'] for i in report.keys()]

    def identify_hirshfield(self):
        # TODO use multwfn to generate hirshfield finger print, then analysze it based on
        # 10.1039/c1ce05763d
        pass
