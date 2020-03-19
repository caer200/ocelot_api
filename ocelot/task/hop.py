import json
import os

import numpy as np
from joblib import Parallel
from joblib import delayed

from ocelot.schema.conformer import DimerCollection
from ocelot.task.zindo import ZindoJob
from ocelot.task.zindo import conver2zindo

"""
for a certain crystal config, get perocalation pathway
"""


class IJK:
    def __init__(self, ijk, klen):
        self.ijk = ijk
        self.klen = klen

    def __eq__(self, other):
        if self.klen != other.klen:
            return False
        i1, j1, k1 = self.ijk
        i2, j2, k2 = other.ijk
        if i1 == j2 and j1 == i2 and k2 == self.klen - k1 - 1:
            return True
        return False


class Hop:

    def __init__(self, config, zindobin=None, zindoctbin=None, zindolib=None, wdir=os.getcwd(), gaussbin=None):
        self.workdir = wdir  # make sure this is an absolute path
        self.zindobin = zindobin
        self.zindoctbin = zindoctbin
        self.zindolib = zindolib
        self.gaussbin = gaussbin
        self.config = config
        self.dimer_array, self.tranv_fcs = self.config.get_dimers_array(maxfold=2, fast=True, symm=False)

    @staticmethod
    def unique_ijk(leni, lenj, lenk):
        """
        given dimer_array from config.get_dimers_array with symm=False
        get unique ijk combinations
        using symmetry dimer_array[i][j][k1] == dimer_array[j][i][k2] where transv_fcs[k1] == -transv_fcs[k2]
        note there is another symmetry in transv_fcs from cart product, transv_fcs[i] == -transv_fcs[len-i-1]
        that is k2 == len-k1-1

        :param leni:
        :param lenj:
        :param lenk:
        :return:
        """
        i_s = list(range(leni))
        j_s = list(range(lenj))
        k_s = list(range(lenk))
        ijks = np.array(np.meshgrid(i_s, j_s, k_s)).T.reshape(-1, 3)
        unique_ijks = []
        for ijk in ijks:
            ijkobj = IJK(ijk, lenk)
            if ijkobj not in unique_ijks:
                unique_ijks.append(ijkobj)
        return [ijkob.ijk for ijkob in unique_ijks]

    @staticmethod
    def screen_dimers(dimer_array, close=True, unique=True):
        """

        :param dimer_array: ijk array as in config.get_dimers_array()
        :param close:
        :param unique: apply symmetry in dimer_array
        :return a list of dimers that are unique and close
        """
        results = []
        if unique:
            ni, nj, nk = dimer_array.shape
            unique_ijks = Hop.unique_ijk(ni, nj, nk)
            for ijk in unique_ijks:
                i, j, k = ijk
                dimer = dimer_array[i, j, k]
                if close:
                    if dimer.is_bone_close and dimer.is_not_identical:
                        results.append(dimer)
                else:
                    if dimer.is_not_identical:
                        results.append(dimer)

        else:
            for dimer in dimer_array.flatten():
                if close:
                    if dimer.is_bone_close and dimer.is_not_identical:
                        results.append(dimer)
                else:
                    if dimer.is_not_identical:
                        results.append(dimer)
        return results

    @staticmethod
    def group_dimers(dimers):
        """
        subgroup into dimercollection such that all dimers in one collection share the same ref_omol

        :param dimers:
        :return: dimercolls: a dict of dimer lists, dimercolls[x] are dimers with i==x
        """
        dimercolls = {}
        for d in dimers:
            i, j, k = [int(idx) for idx in d.label.split('_')]  # I should prob make ijk an attribute
            if i not in dimercolls.keys():
                dimercolls[i] = [d]
            else:
                dimercolls[i].append(d)
        return dimercolls

    def run(self, alldimers, calculator='zindo', njobs=None):
        """
        you should screen alldimers based on geodist before this run

        :param alldimers:
        :param calculator:
        :param njobs:
        :return:
        """
        print('calculator was set to {}'.format(calculator))
        dimercolls = self.group_dimers(alldimers)
        workload = []
        for i in dimercolls.keys():
            dimers = dimercolls[i]
            dimercoll = DimerCollection(dimers)
            dimercoll.to_xyz('{}/dimercoll_{}.xyz'.format(self.workdir, i))
            for idimer in range(len(dimers)):
                lab = dimers[idimer].label
                subwdir = '{}/dimers/{}'.format(self.workdir, lab)
                if calculator != 'geodist':
                    os.system('mkdir -p {}'.format(subwdir))
                mol_A = conver2zindo(dimers[idimer].conformer_ref.pmgmol)
                mol_D = conver2zindo(dimers[idimer].conformer_var.pmgmol)
                xyzstring = dimers[idimer].to_xyzstring()
                unit_job = (lab, subwdir, mol_A, mol_D, xyzstring)
                workload.append(unit_job)
        # run unit jobs
        print('workload has {} unit jobs'.format(len(workload)))
        if calculator == 'zindo':
            if not isinstance(njobs, int):
                njobs = min(os.cpu_count(), len(workload))
            work_results = Parallel(njobs)(delayed(self.zindo_cal)(uj) for uj in workload)
            # total 16 benzene
            # 1 16.880834579467773
            # 2 14.692869424819946
            # 4 17.32540798187256
            # total 3 tipgebw
            # 1 88.71285700798035
            # 3 34.4242103099823
        elif calculator == 'geodist':
            # since alldimers should already have been screened based on geodist, calculator will just return inf
            work_results = [[uj[0], np.inf, np.inf, uj[4]] for uj in workload]
        else:
            raise NotImplementedError('{} is not a valid calculator!'.format(calculator))
        data = {}
        for result in work_results:
            lab, hhti, llti, xyzstring = result
            data[lab] = [hhti, llti, xyzstring]
        with open('{}/hopdata1.json'.format(self.workdir), 'w') as f:
            json.dump(data, f)
        return data

    def zindo_cal(self, unitjob):
        lab, subwdir, mol_A, mol_D, xyzstring = unitjob
        coupling_data, nmo_a, nmo_d = ZindoJob.dimer_run(lab, subwdir, self.zindobin, self.zindoctbin, self.zindolib,
                                                         mol_A, mol_D)
        hhti = ZindoJob.get_hh_coupling(coupling_data, nmo_a, nmo_d)
        llti = ZindoJob.get_ll_coupling(coupling_data, nmo_a, nmo_d)
        return [lab, hhti, llti, xyzstring]

    #
    # def run_zindo(self, alldimers, workdir):
    #     """
    #     Run zindo for all the dimers provided and calculate electronic coupling.
    #     xzy files are written for each of the dimer in the respective folder. xyz file for the dimers is also written.
    #
    #     :param alldimers:
    #     :param workdir: working directory to run zindo
    #     :return: data[dimer.label] = [hhcoupling, llcoupling, xyzstring]
    #     """
    #     data = {}
    #     dimercolls = self.group_dimers(alldimers)
    #     for i in dimercolls.keys():
    #         dimers = dimercolls[i]
    #         dimercoll = DimerCollection(dimers)
    #         dimercoll.to_xyz('{}/dimercoll_{}.xyz'.format(workdir, i))
    #         for idimer in range(len(dimers)):
    #             lab = dimers[idimer].label
    #             subwdir = '{}/dimers/{}'.format(workdir, lab)
    #             os.system('mkdir -p {}'.format(subwdir))
    #
    #             mol_A = conver2zindo(dimers[idimer].conformer_ref.pmgmol)
    #             mol_D = conver2zindo(dimers[idimer].conformer_var.pmgmol)
    #             coupling_data, nmo_a, nmo_d = ZindoJob.dimer_run(lab,
    #                                                              subwdir,
    #                                                              self.zindobin,
    #                                                              self.zindoctbin,
    #                                                              self.zindolib,
    #                                                              mol_A, mol_D)
    #             hhti = ZindoJob.get_hh_coupling(coupling_data, nmo_a, nmo_d)
    #             llti = ZindoJob.get_ll_coupling(coupling_data, nmo_a, nmo_d)
    #             data[lab] = [hhti, llti, dimers[idimer].to_xyzstring()]
    #             dimers[idimer].to_xyz('{}/{}.xyz'.format(subwdir, lab))
    #     with open('{}/hopdata1.json'.format(workdir), 'w') as f:
    #         json.dump(data, f)
    #     return data
    #
    @staticmethod
    def applysym_to_coupling_data(data, dimer_array, workdir):
        leni, lenj, lenk = dimer_array.shape
        symdata = {}
        for lab in data.keys():
            i, j, k = [int(idx) for idx in lab.split('_')]
            symdata[lab] = data[lab]
            symlab = '_'.join([str(index) for index in [j, i, lenk - k - 1]])
            symlabentry = [data[lab][0], data[lab][1], dimer_array[j, i, lenk - k - 1].to_xyzstring()]
            symdata[symlab] = symlabentry
        with open('{}/symdata2.json'.format(workdir), 'w') as f:
            json.dump(data, f)
        return symdata

    @staticmethod
    def get_intger_mesh(x, y, z):
        xs = list(range(-x, x + 1))
        ys = list(range(-y, y + 1))
        zs = list(range(-z, z + 1))
        return np.array(np.meshgrid(xs, ys, zs), dtype=int).T.reshape(-1, 3)

    @staticmethod
    def supercell_proj(dimer_array, transv_fcs, symdata, super_cell_size, motype, workdir):
        """
        we look at a supercell, in which the coupling is represented as augmented_data[mp_i][mp_j][i][j]
        which is the coupling between supercellmesh[mp_i] ith mol and supercellmesh[mp_j] jth mol

        :param motype:
        :param workdir:
        :param dimer_array:
        :param transv_fcs:
        :param symdata:
        :param super_cell_size: 3x1 int tuple e.g. (6, 6, 6)
        :return: supercellmesh, augmented_data
        """
        if motype == 'hh':
            moid = 0
        elif motype == 'll':
            moid = 1
        else:
            moid = 0

        leni, lenj, lenk = dimer_array.shape
        size_x, size_y, size_z = [int(s) for s in super_cell_size]
        super_cell_mesh = Hop.get_intger_mesh(size_x, size_y, size_z)
        mesh_size = len(super_cell_mesh)
        augmented_data = np.zeros((mesh_size, mesh_size, leni, lenj))
        mat_v_x2y = super_cell_mesh[np.newaxis, :, :] - super_cell_mesh[:, np.newaxis, :]
        for mi in range(mesh_size):
            for mj in range(mesh_size):
                # v_mi2mj = super_cell_mesh[mj] - super_cell_mesh[mi]
                v_mi2mj = mat_v_x2y[mi][mj]
                lookup_v_in_transv_fcs = np.argwhere(np.all(transv_fcs == v_mi2mj, axis=1))
                if len(lookup_v_in_transv_fcs) > 0:
                    k = lookup_v_in_transv_fcs[0][0]
                    for i in range(leni):
                        for j in range(lenj):
                            label = '{}_{}_{}'.format(i, j, k)
                            if label in symdata.keys():
                                augmented_data[mi][mj][i][j] = symdata[label][moid]
                                augmented_data[mj][mi][j][i] = symdata[label][moid]

        with open('{}/supercelldata3.json'.format(workdir), 'w') as f:
            json.dump(augmented_data.tolist(), f)
        return super_cell_mesh, augmented_data

    @staticmethod
    def hopping(supercell_data, supercell_mesh, cutoff, workdir, motype):
        """
        a molecule is defined by meshpoint in the supercell mesh (mp_i) and mol index (mol_i) in a cell
        hopdata[x] gives a list of molecules (represented by [mp_i, mol_i])
        that are connected to the molecule [mp_i=0, mol_i=x], including the molecule itself

        :param motype:
        :param supercell_mesh: 3xn meshpoints
        :param supercell_data:
        :param cutoff: in meV
        :param workdir:
        :return:
        """
        meshsize = len(supercell_mesh)
        number_of_mols = len(supercell_data[0][0])

        def conneted_neighbors(mp_i, mol_i):
            list_of_connected = []
            for mp_j in range(meshsize):
                for mol_j in range(number_of_mols):
                    ti = supercell_data[mp_i][mp_j][mol_i][mol_j]
                    if abs(ti) > cutoff:
                        list_of_connected.append((mp_j, mol_j))
            return list_of_connected

        network = {}
        for mol_init_i in range(number_of_mols):
            # print('start {}'.format(mol_init_i))
            mp_i_init = int((meshsize - 1) / 2)  # s.t. mesh[mp_i_init] = (0, 0, 0)
            in_network = [(mp_i_init, mol_init_i)]
            pointer = 0
            while len(in_network) != pointer:
                focus_mol = in_network[pointer]
                # print('looking at mp {} mol {}'.format(*in_network[pointer]))
                for connected in conneted_neighbors(*focus_mol):
                    if connected not in in_network:
                        # print('find connected {}'.format(connected))
                        in_network.append(connected)
                pointer += 1
            network[mol_init_i] = in_network
        with open('{}/{}_netdata4.json'.format(workdir, motype), 'w') as f:
            json.dump(network, f)
        return network

    def get_hopping_network_s1(self):
        dimers_to_be_used = self.screen_dimers(self.dimer_array, close=True, unique=True)
        # this gives 'hopdata1.json'
        # hopdata = self.run_zindo(dimers_to_be_used, self.workdir)
        hopdata = self.run(dimers_to_be_used, 'zindo')
        # this gives 'symdata2.json'
        symdata = self.applysym_to_coupling_data(hopdata, self.dimer_array, self.workdir)
        # this gives 'supercelldata3.json'
        return hopdata, symdata

    def get_hopping_network_s2(self, symdata, cutoff, supercell=(1, 1, 1), motype='hh'):
        mesh, augdata = self.supercell_proj(self.dimer_array, self.tranv_fcs, symdata, super_cell_size=supercell,
                                            motype=motype, workdir=self.workdir)
        # this gives 'll_netdata4.json' or 'hh_netdata4.json'
        network = self.hopping(augdata, mesh, cutoff, self.workdir, motype)
        return mesh, augdata, network

    # # this is an old implementaion with tuple as key for augdata, very slow
    # augmented_data = {}
    # for meshpoint_x in super_cell_mesh:
    #     augmented_data[meshpoint_x] = {}
    #     for meshpoint_y in super_cell_mesh:
    #         augmented_data[meshpoint_x][meshpoint_y] = {}
    #         v_x2y = tuple(meshpoint_y[i] - meshpoint_x[i] for i in range(3))
    #         try:
    #             k = transv_fcs.index(v_x2y)
    #         except ValueError:
    #             k = None
    #         if k is not None:
    #             for i in range(leni):
    #                 augmented_data[meshpoint_x][meshpoint_y][i] = {}
    #                 for j in range(lenj):
    #                     label = '{}_{}_{}'.format(i, j, k)
    #                     if label in symdata.keys():
    #                         augmented_data[meshpoint_x][meshpoint_y][i][j] = symdata[label]
    #                     else:
    #                         augmented_data[meshpoint_x][meshpoint_y][i][
    #                             j] = 0.0  # self coupling is already included here
    #         else:
    #             for i in range(leni):
    #                 augmented_data[meshpoint_x][meshpoint_y][i] = {}
    #                 for j in range(lenj):
    #                     augmented_data[meshpoint_x][meshpoint_y][i][j] = 0.0
    # return augmented_data

    # @staticmethod
    # def hopping(supercell_data, cutoff, workdir):
    #     """
    #
    #
    #     :param supercell_data:
    #     :param cutoff: in meV
    #     :param workdir:
    #     :return:
    #     """
    #
    #     mps = list(supercell_data.keys())
    #     js = list(supercell_data[mps[0]][mps[0]][0].keys())
    #
    #     def conneted_neighbors(mpi):
    #         mp_x = mpi[:3]
    #         i = mpi[3]
    #         list_of_mp_y_j = []
    #         for mp_y in mps:
    #             for j in js:
    #                 ti = supercell_data[mp_x][mp_y][i][j]
    #                 if abs(ti) > cutoff:
    #                     list_of_mp_y_j.append((*mp_y, j))
    #         return list_of_mp_y_j
    #
    #     network = {}
    #     for init_i in js:
    #         mpi_init = (0, 0, 0, init_i)
    #         in_network = [mpi_init]
    #         pointer = 0
    #         while len(in_network) != pointer:
    #             focus_mpi = in_network[pointer]
    #             for connected in conneted_neighbors(focus_mpi):
    #                 if connected not in in_network:
    #                     in_network.append(connected)
    #             pointer += 1
    #         network[init_i] = in_network
    #     with open('{}/netdata4.json'.format(workdir), 'w') as f:
    #         json.dump(network, f)
    #     return network
