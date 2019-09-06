import numpy as np
import random
import sys
from copy import deepcopy
from api.schema.backbone import Backbone
from api.schema.Backbone import Backbone
from api.schema.msitelist import Sitelist
from scipy.spatial.distance import cdist
from api.routines.geometry import rotate_along_axis, unify, rotation_matrix, norm, angle_btw


class Boxfitter:

    def __init__(self, mol, pbc):
        self.mol = deepcopy(mol)
        self.pbc = np.array(pbc)
        a, b, c = self.pbc
        self.transvs = [a, b, c, a - b, a - c, b - c, a + b, a + c, b + c]

    def tooclose(self, config, cutoff=3):
        coordmat = config.get_coordmat()
        for v in self.transvs:
            distmat = cdist(coordmat, coordmat + v, 'euclidean')
            if np.min(distmat) < cutoff:
                return True
        return False

    def gen_bone_configs(self, steps=180):
        backbone = Backbone.orient(self.mol.backbone)  # deepcopy, geoc=origin vp=x vq=y vo=z
        bone_configs = []
        step_size = 180.0 / steps
        nconfigs = (steps * 2 + 1) ** 2 * (steps + 1)
        config_coords = np.zeros((nconfigs, len(backbone.sites), 3))
        config_pqos = np.zeros((nconfigs, 3, 3))
        iconfig = 0
        rotmat_z = rotation_matrix(axis='z', theta=step_size, thetaunit='degree')
        vp = np.array([1, 0, 0])
        vq = np.array([0, 1, 0])
        for i in range(steps * 2 + 1):
            if i != 0:
                backbone.rotate_with_matrix(rotmat_z, (0, 0, 0))
                vp = rotate_along_axis(vp, axis=(0, 0, 1), theta=step_size, thetaunit='degree')
                vq = rotate_along_axis(vq, axis=(0, 0, 1), theta=step_size, thetaunit='degree')
            rotmat_vq = rotation_matrix(axis=vq, theta=step_size, thetaunit='degree')
            rotmat_vq_180 = rotation_matrix(axis=vq, theta=180, thetaunit='degree')
            for j in range(steps + 1):
                if j != 0:
                    backbone.rotate_with_matrix(rotmat_vq, (0, 0, 0))
                    vp = rotate_along_axis(vp, axis=vq, theta=step_size, thetaunit='degree')
                rotmat_vp = rotation_matrix(axis=vp, theta=step_size, thetaunit='degree')
                for k in range(steps * 2 + 1):
                    if k != 0:
                        backbone.rotate_with_matrix(rotmat_vp, (0, 0, 0))
                        vq = rotate_along_axis(vq, axis=vp, theta=step_size, thetaunit='degree')
                    vo = np.cross(vp, vq)
                    config_coords[iconfig] = backbone.get_coordmat()
                    config_pqos[iconfig] = [vp, vq, vo]
                    iconfig += 1
            backbone.rotate_with_matrix(rotmat_vq_180, (0, 0, 0))
        names = [s.element.name for s in backbone.sites]
        ids = [s.id for s in backbone.sites]
        for coords in config_coords:
            bone_configs.append(Sitelist.from_coordmat(coords, names, ids))
        return bone_configs, config_pqos

    def addsj(self, config, sidechains, bone_vpqo):
        for sc in sidechains:
            bj = config.get_site(sc.bone_joint.id)
            v = sc.side_joint.coords - sc.bone_joint.coords
            vcoord = np.matmul(v, self.mol.backbone.pqomat.T)
            sjs = deepcopy(sc.side_joint)
            sjs.coords = np.matmul(vcoord, bone_vpqo) + bj.coords
            config.sites.append(sjs)
        return config

    def addtosite(self, config, steps=180):
        stepsize = 180.0 / steps
        names = [s.element.name for s in config.sites]
        ids = [s.id for s in config.sites]
        # allids = list(range(len(self.mol.msites)))
        configids = [s.id for s in config]
        current_site = None
        for s in config:
            s_origin = self.mol.get_site(s.id)
            if hasattr(s_origin, 'higher_rank_nbs'):
                if len(s_origin.higher_rank_nbs) > 0:
                    if s_origin.higher_rank_nbs[0].id not in configids:
                        current_site = s
                        break
        if current_site is None:
            return [config]
        else:
            configs = []
            current_site_in_mol = self.mol.get_site(current_site.id)
            root_site = config.get_site(current_site_in_mol.lower_rank_nbs[0].id)
            v_root_to_current = unify(current_site.coords - root_site.coords)
            # TODO right now only sp3c-like elements are supported
            if current_site_in_mol.hybrid == 'sp3':
                configs_coords = np.zeros((steps * 2, len(configids) + 3, 3))
                leaf1, leaf2, leaf3 = current_site_in_mol.higher_rank_nbs
                names += [leaf1.element.name, leaf2.element.name, leaf3.element.name]
                ids += [leaf1.id, leaf2.id, leaf3.id]
                b1 = norm(leaf1.coords - current_site_in_mol.coords)
                b2 = norm(leaf2.coords - current_site_in_mol.coords)
                b3 = norm(leaf3.coords - current_site_in_mol.coords)
                if abs(angle_btw(np.array([1, 0, 0]), v_root_to_current)) < 1e-3 or abs(
                        angle_btw(np.array([1, 0, 0]), v_root_to_current) - np.pi) < 1e-3:
                    vertical_to_vrc = np.cross((1, 0, 0), v_root_to_current)
                else:
                    vertical_to_vrc = np.cross((0, 1, 0), v_root_to_current)
                leaf1_coords = rotate_along_axis(-v_root_to_current, vertical_to_vrc, 109.5,
                                                 'degree') * b1 + current_site.coords
                leaf2_coords = rotate_along_axis(unify(leaf1_coords - current_site.coords), v_root_to_current, 120,
                                                 'degree') * b2 + current_site.coords
                leaf3_coords = rotate_along_axis(unify(leaf1_coords - current_site.coords), v_root_to_current, 240,
                                                 'degree') * b3 + current_site.coords
                for i in range(steps * 2):
                    if i != 0:
                        leaf1_coords = rotate_along_axis(leaf1_coords - current_site.coords, v_root_to_current,
                                                         stepsize, 'degree') + current_site.coords
                        leaf2_coords = rotate_along_axis(leaf2_coords - current_site.coords, v_root_to_current,
                                                         stepsize, 'degree') + current_site.coords
                        leaf3_coords = rotate_along_axis(leaf3_coords - current_site.coords, v_root_to_current,
                                                         stepsize, 'degree') + current_site.coords
                    configs_coords[i] = np.vstack((config.get_coordmat(), leaf1_coords, leaf2_coords, leaf3_coords))
            elif current_site_in_mol.hybrid == 'sp':
                configs_coords = np.zeros((1, len(configids) + 1, 3))
                leaf1 = current_site_in_mol.higher_rank_nbs[0]
                names += [leaf1.element.name]
                ids += [leaf1.id]
                b1 = norm(leaf1.coords - current_site_in_mol.coords)
                leaf1_coords = v_root_to_current * b1 + current_site.coords
                configs_coords[0] = np.vstack((config.get_coordmat(), leaf1_coords))
            else:
                sys.exit('undefined hybrid!')
            for coords in configs_coords:
                configs.append(Sitelist.from_coordmat(coords, names, ids))
            return configs

    def grow(self, bone_config, bone_pqo, steps, expansion, itag):
        c = self.survive([bone_config])
        if len(c) == 0:
            return c
        c = [self.addsj(bone_config, self.mol.sidechains, bone_pqo)]
        for i in range(expansion):
            if len(c) == 0:
                return c
            r = []
            for config in c:
                r += self.addtosite(config, steps)
            c = self.survive(r)
            print('expansion', i, 'sur', len(c))
            if len(c) > 0 and i > 10:
                j = 0
                for con in c:
                    con.to_poscar('c-{}-{}-{}.poscar'.format(itag, j, i), self.pbc)
                    j += 1
        return c

    def survive(self, configs):
        sur = []
        for c in configs:
            if not self.tooclose(c):
                sur.append(c)
        if len(sur) > 500:
            sur = random.sample(sur, 500)
        return sur
