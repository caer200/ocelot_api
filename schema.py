from scipy.spatial.distance import pdist, squareform
import pyny3d.geoms as pyny
import json
import subprocess
import itertools
from collections import Counter
from copy import deepcopy
from geo import *


class Site:
    def __init__(self, ele, coords):
        self.ele = ele
        self.coords = np.array(coords)

    def dist(self, other):
        return np.linalg.norm(self.coords - other.coords)

    @property
    def radii(self):
        return covalent_radii[self.ele]

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    def as_dict(self):
        d = {}
        d['ele'] = self.ele
        d['coords'] = self.coords.tolist()  # np array not json serializable
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d['ele'], d['coords'])


def coordmat(sites):
    return np.array([s.coords for s in sites])


def triangle_area(sites):
    if len(sites) == 3:
        triangle = pyny.Polygon(coordmat(sites))
        return triangle.get_area()
    else:
        sys.exit('build a triangle with nsites != 3')


def get_nelec(sites):
    nelec = 0
    for s in sites:
        nelec += atomic_numbers[s.ele]
    return nelec


def get_nbmap(sites, cutoff_ratio=1.5, absolute=1.99, cutoff_type='absolute'):
    nsites = len(sites)
    distmat = squareform(pdist(coordmat(sites)))
    bmat = np.ones((nsites, nsites), dtype=bool)
    for i in range(nsites):
        bmat[i][i] = False
        for j in range(i + 1, nsites):
            if cutoff_type == 'absolute':
                cutoff = absolute
            else:
                isite = sites[i]
                jsite = sites[j]
                cutoff = (isite.radii + jsite.radii) * cutoff_ratio
            if distmat[i][j] > cutoff:
                bmat[i][j] = False
                bmat[j][i] = False
    nbmap = {}
    for i in range(nsites):
        nbs = []
        for j in range(nsites):
            if bmat[i][j] and j != i:
                nbs.append(j)
        nbmap[i] = nbs
    return nbmap


class Defect:
    def __init__(self, sites, codename=None):
        self.sites = [s for s in sites if s.ele != 'V']
        self.codename = str(codename).split('/')[-1]
        self.nbmap = get_nbmap(self.sites)

    def __eq__(self, other):
        if len(self.sites) != len(other.sites):
            return False

        if self.stoi != other.stoi:
            return False

        if self.edge_stoi != other.edge_stoi:
            return False

        # if len(self.hetero_sites) == 2:
        #     h1, h2 = self.hetero_sites
        #     oh1, oh2 = other.hetero_sites
        #     if abs(h1.dist(h2) - oh1.dist(oh2)) > 1e-3:
        #         return False
        #
        # if len(self.edge_sites) == 2:
        #     h1, h2 = self.edge_sites
        #     oh1, oh2 = other.edge_sites
        #     if abs(h1.dist(h2) - oh1.dist(oh2)) > 1e-3:
        #         return False
        #
        # if len(self.hetero_sites) == 3:
        #     area_self = triangle_area(self.hetero_sites)
        #     area_other = triangle_area(other.hetero_sites)
        #     if abs(area_other - area_self) > 1e-3:
        #         return False
        #
        # if len(self.edge_sites) == 3:
        #     area_self = triangle_area(self.edge_sites)
        #     area_other = triangle_area(other.edge_sites)
        #     if abs(area_other - area_self) > 1e-3:
        #         return False

        if self.doistring != other.doistring:
            return False

        return True

    @property
    def nelec(self):
        return get_nelec(self.sites)

    @classmethod
    def addh(cls, d):
        ids_tobeadded = d.unsaturated_carbon_ids
        hsites = []
        for i in ids_tobeadded:
            s_tobeadded = d.sites[i]
            nbs_id = d.nbmap[i]
            s1, s2 = [d.sites[j] for j in nbs_id]
            v1 = np.array(s1.coords) - np.array(s_tobeadded.coords)
            v2 = np.array(s2.coords) - np.array(s_tobeadded.coords)
            nv = -(v1 + v2)
            nv = nv / np.linalg.norm(nv)
            hsite_coords = s_tobeadded.coords + nv
            hsite = Site('H', hsite_coords.tolist())
            hsites.append(hsite)
        nsites = deepcopy(d.sites)
        nsites += hsites
        return cls(nsites, d.codename + 'H')

    def addsorption_coords(self, extend):
        dopant_ids = [sid for sid in range(len(self.sites)) if self.sites[sid].ele in ('S', 'N')]
        centers = []  # TODO may be useful
        # if len(dopant_ids) >= 2:
        #     combs = itertools.combinations(dopant_ids, 2)
        #     for c in combs:
        #         dist = np.linalg.norm(np.array(self.sites[c[0]].coords) - np.array(self.sites[c[1]].coords))
        #         if dist < 3.5:
        #             center_coords = (np.array(self.sites[c[0]].coords) + np.array(self.sites[c[1]].coords)) / 2
        #             centers.append(center_coords)

        # ads_ids = list(set(dopant_ids))

        dopant_nbs = []
        for i in dopant_ids:
            dopant_nbs += self.nbmap[i]
        ads_ids = list(set(dopant_ids + dopant_nbs))

        ads_ids_noh = []
        for i in ads_ids:
            if not self.sites[i].ele == 'H':
                ads_ids_noh.append(i)
        ads_ids = ads_ids_noh

        # add chem env criteria
        ads_ids_chemenvdiff = []
        chemenvs = []
        for i in ads_ids:
            chemenv = self.chem_env(i, self.nbmap, max_level=6)
            if chemenv not in chemenvs:
                chemenvs.append(chemenv)
                ads_ids_chemenvdiff.append(i)
        ads_ids = ads_ids_chemenvdiff

        ads_coords = [self.sites[i].coords for i in ads_ids] + centers
        ads_coords = np.array(ads_coords)
        ads_coords = [extend * coords / np.linalg.norm(coords) + coords for coords in ads_coords]
        return ads_coords

    @property
    def doistring(self):
        return '@'.join(sorted([self.chem_env(siteid, self.nbmap, max_level=4) for siteid in range(len(self.sites))]))

    @property
    def hetero_or_edge_ids(self):
        idlst = [sid for sid in range(len(self.sites)) if self.sites[sid].ele != 'C' and self.sites[sid].ele != 'H']
        idlst += self.edge_ids
        return idlst

    @property
    def unsaturated_carbon_ids(self):
        return [i for i in self.edge_ids if self.sites[i].ele == 'C']

    def chem_env(self, siteid, nbmap, max_level=2):  # percolation 2 times
        nbs_per_percolation = {}  # nbs_per_percolation[1] = an id list of nbs, nbs_per_percolation[2] = an id list of 2nd nbs
        nbs_per_percolation[0] = [siteid]
        for i in range(max_level):
            nblist = []
            included = []
            for j in range(i + 1):
                included += nbs_per_percolation[j]
            strating_ids = nbs_per_percolation[i]
            for si in strating_ids:
                nblist += [nbid for nbid in nbmap[si] if nbid not in included]
            nbs_per_percolation[i + 1] = nblist

        doi = '{}{}_'.format(0, self.sites[siteid].ele)
        for k in range(max_level + 1):
            elements = sorted([self.sites[nb].ele for nb in nbs_per_percolation[k]])
            doi += '-'
            doi += ''.join(elements)
        return doi

    @property
    def hetero_sites(self):
        return [s for s in self.sites if s.ele != 'C' and s.ele != 'H']

    @property
    def hydrogen_sites(self):
        return [s for s in self.sites if s.ele == 'H']

    @property
    def hetero_distances(self):
        dists = []
        for ij in itertools.combinations(range(len(self.hetero_sites)), 2):
            i, j = ij
            dists.append(self.hetero_sites[i].dist(self.hetero_sites[j]))
        if len(dists) == 0:
            dists = [0]
        return dists

    @property
    def edge_sites(self):
        return [self.sites[i] for i in self.edge_ids]

    @property
    def edge_stoi(self):
        c = Counter([s.ele for s in self.edge_sites])
        uniques = sorted(c.keys())
        formula = ''
        for u in uniques:
            formula += '{}{} '.format(u, c[u])
        return formula

    @property
    def edge_ids(self):  # 2-coord C/S/N or C in a C-H bond
        edge = []
        nsites = len(self.sites)
        for i in range(nsites):
            if len(self.nbmap[i]) == 2 and self.sites[i].ele in ('C', 'S', 'N'):
                edge.append(i)
            elif len(self.nbmap[i]) == 3 and self.sites[i].ele == 'C' and 'H' in [self.sites[ins].ele for ins in
                                                                                  self.nbmap[i]]:
                edge.append(i)
        return edge

    @property
    def stoi(self):
        c = Counter([s.ele for s in self.sites])
        uniques = sorted(c.keys())
        formula = ''
        for u in uniques:
            formula += '{}{} '.format(u, c[u])
        return formula

    def as_dict(self):
        def_dict = dict()
        def_dict['sites'] = siteobjs2sitedicts(self.sites)
        def_dict['codename'] = self.codename
        return def_dict

    @classmethod
    def from_dict(cls, def_dict):
        sites = sitedicts2siteobjs(def_dict['sites'])
        codename = def_dict['codename']
        return cls(sites, codename)

    @classmethod
    def from_file(cls, fn, codename, caltype='rlx'):
        names = fn.split('/')[-1].split('.')
        if len(names) < 2:
            sys.exit('failed to creat defect obj from {} as no extension found'.format(fn))
        extension = names[-1]
        if extension == 'mop':
            sites = mop2siteobjs(fn)
            return cls(sites, codename=codename)
        elif extension == 'xyz':
            sites = xyz2siteobjs(fn)
            return cls(sites, codename=codename)
        elif extension == 'out':
            d = out2out_dict(fn, codename, caltype=caltype)
            sitedicts = d['sites']
            sites = sitedicts2siteobjs(sitedicts)
            return cls(sites, codename=codename)
        else:
            sys.exit('failed to creat defect obj from {} as extension not recognized'.format(fn))

    def to_string(self, extension, mopheadline=None):
        return siteobjs2string(self.sites, extension, self.codename, mopheadline)

    def to(self, fn, extension, mopheadline=None):
        siteobjs2file(self.sites, fn, extension, self.codename, mopheadline)

    @property
    def nN(self):
        return len([s for s in self.sites if s.ele == 'N'])

    @property
    def nH(self):
        return len([s for s in self.sites if s.ele == 'H'])

    @property
    def nS(self):
        return len([s for s in self.sites if s.ele == 'S'])

    @property
    def legit(self):
        """
        1. for any atom, # of neighbors >= 1
        2. for sulfur atoms, # of neighbors == 2
        3. for carbon atoms, # of neighbors == 2, 3
        4. for nitrogen atoms, # of neighbors == 2, 3
        5. no S-S, N-S, N-N bonds
        6. \# of N + # of S should be less than 4
        :return:
        """
        nsites = len(self.sites)
        nN = len([s for s in self.sites if s.ele == 'N'])
        nS = len([s for s in self.sites if s.ele == 'S'])
        if nN + nS >= 4:
            self.nolegit='vio6'
            return False  # violating 6.
        for i in range(nsites):
            s = self.sites[i]
            nbs = self.nbmap[i]
            if len(nbs) < 1:
                self.nolegit='vio1'
                return False  # violating 1.
            if s.ele == 'S':
                if len(nbs) != 2:
                    self.nolegit='vio2'
                    return False  # violating 2.
                if any([self.sites[j].ele in ('N', 'S') for j in nbs]):
                    self.nolegit='vio5'
                    return False  # violating 5.
            if s.ele == 'H' and len(nbs) != 1:
                self.nolegit='vioh'
                return False
            if s.ele == 'C' and len(nbs) not in (2, 3):
                self.nolegit='vio3'
                return False  # violating 3.
            if s.ele == 'N' and len(nbs) not in (2, 3):
                self.nolegit='vio4'
                return False  # violating 4.
            if s.ele in ['N', 'S'] and any([self.sites[j].ele in ('N', 'S') for j in nbs]):
                self.nolegit='vio5'
                return False  # violating 5.
        return True

    @classmethod
    def substitution(cls, site_ids, dopants, basesites, codename=None):
        for i in range(len(site_ids)):
            dopant = dopants[i]
            siteid = site_ids[i]
            basesites[siteid].ele = dopant
        return cls(basesites, codename)

    @staticmethod
    def dumplist(defect_list, fn):
        json_list = [d.as_dict() for d in defect_list]
        with open(fn, 'w') as f:
            json.dump(json_list, f)

    @staticmethod
    def loadlist(fn):
        with open(fn, 'r') as f:
            defectlist = json.load(f)
        for i in range(len(defectlist)):
            defectlist[i] = Defect.from_dict(defectlist[i])
        return defectlist

    def to_mopin(self, mopwd, mopheadline, caltype):
        job_codename = self.codename + '-' + caltype
        mopin = '{}/{}.mop'.format(mopwd, job_codename)
        mopinstring = siteobjs2file(self.sites, fn=mopin, extension='mop', codename=job_codename,
                                    mopheadline=mopheadline)
        return mopinstring, mopin, job_codename

    def cal_thermal(self, mopwd, bincmd, mopheadline):
        out_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, self.codename, mopheadline,
                                                          caltype='thermal')
        return out_dict, job_codename, mopin, mopout

    # def relax(self, mopwd, bincmd, codename, mopheadline_rlx):
    #     arc_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, codename, mopheadline_rlx,
    #                                                       suffix='rlx')  # tot, hof, sites, fstring
    #     arc_dict['codename'] = job_codename
    #     relaxed_sites = sitedicts2siteobjs(arc_dict['sites'])
    #     relaxed_defect = Defect(relaxed_sites, job_codename)
    #     return arc_dict, mopin, mopout, relaxed_defect
    #
    # def cal(self, mopwd, bincmd, codename, mopheadline_rlx, mopheadline_thermal):
    #     arc_dict0, mopin0, mopout0, relaxed_defect0 = self.relax(mopwd, bincmd, codename, mopheadline_rlx,)
    #     arc_dict1, mopin1, mopout1, relaxed_defect1 = relaxed_defect0.relax(mopwd, bincmd, codename, mopheadline_rlx)
    #     arc_dict2, mopin2, mopout2, relaxed_defect2 = relaxed_defect1.relax(mopwd, bincmd, codename, mopheadline_rlx)
    #     arc_dict, job_codename, mopin, mopout = mopac_run(relaxed_defect2.sites, mopwd, bincmd, codename,
    #                                                       mopheadline_thermal, suffix='thermal')
    #     fin_defect = relaxed_defect2
    #     fin_defect.codename = codename + '-relaxed'
    #     gnorm = arc_dict['gnorm']
    #     tot = arc_dict2['tot']
    #     hof = arc_dict['hof']
    #     entropy = arc_dict['entropy']
    #     return tot, hof, entropy, gnorm, fin_defect


class Adsorbate:

    def __init__(self, sites, codename, tot, hof, ads_id):
        self.sites = sites
        self.speciesname = codename
        self.codename = '{}*{}'.format(codename, ads_id)
        self.tot = tot
        self.hof = hof
        self.ads_id = int(ads_id)

    @property
    def ads_site(self):
        return self.sites[int(self.ads_id)]

    def zero(self):
        ref_site_coords = deepcopy(self.ads_site.coords)
        for i in range(len(self.sites)):
            self.sites[i].coords = self.sites[i].coords - ref_site_coords

    @property
    def axis(self):
        end1 = self.ads_site
        if 'cooh*' in self.codename or 'cooneg*' in self.codename:
            end2 = self.sites[1].coords + self.sites[2].coords
        elif 'chooh*' in self.codename:
            end2 = self.sites[1].coords
        elif 'co*' in self.codename:
            end2 = self.sites[1].coords
        else:
            end2 = sorted(self.sites, key=lambda x: np.linalg.norm(x.coords - self.ads_site.coords))[-1]
            end2 = end2.coords
        return end2 - end1.coords

    def as_dict(self):
        d = {}
        d['sites'] = siteobjs2sitedicts(self.sites)
        d['codename'] = self.codename
        d['tot'] = self.tot
        d['hof'] = self.hof
        return d

    @classmethod
    def from_dict(cls, d):
        sites = sitedicts2siteobjs(d['sites'])
        codename, ads_id = d['codename'].split('*')[:2]
        return cls(sites, codename, d['tot'], d['hof'], ads_id)

    @classmethod
    def from_file(cls, fn, codename, ads_id, suffix='rlx'):
        names = fn.split('/')[-1].split('.')
        if len(names) < 2:
            sys.exit('failed to creat defect obj from {} as no extension found'.format(fn))
        extension = names[-1]
        if extension == 'mop':
            sites = mop2siteobjs(fn)
            return cls(sites, codename, ads_id=ads_id, tot=None, hof=None)
        elif extension == 'xyz':
            sites = xyz2siteobjs(fn)
            return cls(sites, codename, ads_id=ads_id, tot=None, hof=None)
        elif extension == 'out':
            d = out2out_dict(fn, codename, caltype=suffix)
            sites = sitedicts2siteobjs(d['sites'])
            hof = d['hof']
            tot = d['tot']
            return cls(sites, codename, tot, hof, ads_id)
        else:
            sys.exit('failed to creat adsorbate obj from {} as extension not recognized'.format(fn))

    @classmethod
    def place(cls, adsobj, vref, origin, angle_ref, angle_perpen):
        if 'chooh' not in adsobj.codename:
            ads = deepcopy(adsobj)
            ads.zero()
            theta = angle_btw(vref, ads.axis, output='radian')
            v_perpen = unify(np.cross(vref, ads.axis))
            for i in range(len(ads.sites)):
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, v_perpen, -theta, thetaunit='radian')
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, v_perpen, angle_perpen, thetaunit='degree')
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, vref, angle_ref, thetaunit='degree')
            ads_site_coords = np.array(ads.sites[adsobj.ads_id].coords)
            v_trans = np.array(origin) - ads_site_coords
            for i in range(len(ads.sites)):
                ads.sites[i].coords = ads.sites[i].coords + v_trans
            return ads
        else:
            ads = deepcopy(adsobj)
            ads.zero()
            theta = angle_btw(vref, ads.axis, output='radian')
            v_perpen = unify(np.cross(vref, ads.axis))
            for i in range(len(ads.sites)):
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, v_perpen, -theta, thetaunit='radian')
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, v_perpen, angle_perpen, thetaunit='degree')
                ads.sites[i].coords = rotate_along_axis(ads.sites[i].coords, vref, angle_ref, thetaunit='degree')
            ads_site_coords = np.array(ads.sites[adsobj.ads_id].coords)
            v_trans = np.array(origin) - ads_site_coords
            for i in range(len(ads.sites)):
                ads.sites[i].coords = ads.sites[i].coords + v_trans

            c0 = ads.sites[0].coords
            o3 = ads.sites[3].coords
            o1 = ads.sites[1].coords
            v_co3 = o3 - c0
            v_co1 = o1 - c0
            angle_1 = np.pi / 2.0 - angle_btw(v_co3, c0, 'radian')
            # angle_1 = 0
            for i in range(len(ads.sites)):
                vci = ads.sites[i].coords - c0
                nvci = rotate_along_axis(vci, v_co1, angle_1, thetaunit='radian')
                ads.sites[i].coords = c0 + nvci
            return ads

    @property
    def legit(self):
        nsites = len(self.sites)
        nbmap = get_nbmap(self.sites)
        for i in range(nsites):
            nbs = nbmap[i]
            if len(nbs) < 1:
                return False  # violating 1.
        visited = []
        block_list = []
        indices = range(nsites)
        while len(visited) != nsites:
            unvisited = [idx for idx in indices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            pointer = 0
            while pointer != len(block):
                block += [idx for idx in nbmap[block[pointer]] if idx not in block and idx not in visited]
                visited.append(block[pointer])
                pointer += 1
            block_list.append(block)
        if len(block_list) != 1:
            return False
        return True

    # def relax(self, mopwd, bincmd, codename, mopheadline):
    #     arc_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, codename, mopheadline,
    #                                                       suffix='rlx')
    #     tot = arc_dict['tot']
    #     hof = arc_dict['hof']
    #     relaxed_sites = sitedicts2siteobjs(arc_dict['sites'])
    #     relaxed_ads = Adsorbate(relaxed_sites, job_codename, tot, hof, self.ads_id)
    #     return tot, relaxed_ads

    # def cal(self, mopwd, bincmd, codename, mopheadline_rlx, mopheadline_thermal):
    #     arc_dict_rlx, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, codename, mopheadline_rlx,
    #                                                       suffix='rlx')
    #     tot = arc_dict_rlx['tot']
    #     hof = arc_dict_rlx['hof']
    #     relaxed_sites = sitedicts2siteobjs(arc_dict_rlx['sites'])
    #     relaxed_ads = Adsorbate(relaxed_sites, job_codename, tot, hof, self.ads_id)
    #     arc_dict_thermal, job_codename, mopin, mopout = mopac_run(relaxed_ads.sites, mopwd, bincmd, codename,
    #                                                       mopheadline_thermal, suffix='thermal')
    #     hof = arc_dict_thermal['hof']
    #     entropy = arc_dict_thermal['entropy']
    #     return tot, hof, entropy, arc_dict_rlx, arc_dict_thermal,


class Config:
    def __init__(self, defect, adsorbate, defect_adsid, angle_ref, angle_perpen, extend=1.8, apply_placing=True):
        self.codename = '{}-{}-{}-{}-{}-{}'.format(defect.codename, adsorbate.codename,
                                                   defect_adsid, angle_ref, angle_perpen, extend)
        self.defect_adsid = defect_adsid
        self.angle_ref = angle_ref
        self.angle_perpen = angle_perpen
        self.extend = extend
        self.defect = deepcopy(defect)
        if apply_placing:
            vref = unify(defect.addsorption_coords(extend=extend)[defect_adsid])
            self.adsorbate = Adsorbate.place(
                adsobj=adsorbate,
                vref=vref,
                origin=self.defect.addsorption_coords(extend=extend)[defect_adsid],
                angle_ref=angle_ref,
                angle_perpen=angle_perpen
            )
        else:
            self.adsorbate = deepcopy(adsorbate)
        self.sites = self.adsorbate.sites + self.defect.sites  # always ads sites come first

    def as_dict(self):
        d = {}
        d['codename'] = self.codename
        d['defect'] = self.defect.as_dict()
        d['adsorbate'] = self.adsorbate.as_dict()
        d['defect_adsid'] = self.defect_adsid
        d['angle_ref'] = self.angle_ref
        d['angle_perpen'] = self.angle_perpen
        return d

    @classmethod
    def from_dict(cls, d, apply_placing=True):
        defect = Defect.from_dict(d['defect'])
        adsorbate = Adsorbate.from_dict(d['adsorbate'])
        defect_adsid = d['defect_adsid']
        angle_ref = d['angle_ref']
        angle_perpen = d['angle_perpen']
        return cls(defect, adsorbate, defect_adsid, angle_ref, angle_perpen, apply_placing=apply_placing)

    @classmethod
    def from_outfile(cls, outfile, def_codename, adsorbate_mol, defect_adsid, angle_ref, angle_perpen, caltype,
                     inherit=False):
        out_dict = out2out_dict(outfile, None, caltype=caltype)
        sites = sitedicts2siteobjs(out_dict['sites'])
        natoms_ads = len(adsorbate_mol.sites)
        ads_sites = sites[:natoms_ads]
        defect_sites = sites[natoms_ads:]
        defect = Defect(defect_sites, def_codename)
        # print(len(defect.addsorption_coords()), 'defect bug_fromfile', defect.addsorption_coords())
        #
        # debugsites = defect.sites + [Site('La', coords) for coords in defect.addsorption_coords()]
        # siteobjs2file(debugsites, 'bug3la.xyz', extension='xyz', codename='bug1')
        # print(len(defect.addsorption_coords()), 'bug3', )

        adsorbate = Adsorbate(ads_sites, adsorbate_mol.codename, adsorbate_mol.tot, adsorbate_mol.hof,
                              adsorbate_mol.ads_id)
        return cls(defect, adsorbate, defect_adsid, angle_ref, angle_perpen, apply_placing=False)

    def to_string(self, extension, mopheadline=None):
        return siteobjs2string(self.sites, extension, self.codename, mopheadline)

    def to(self, fn, extension, mopheadline=None):
        siteobjs2file(self.sites, fn, extension, self.codename, mopheadline)

    @property
    def legit(self):
        nsites = len(self.sites)
        nbmap = get_nbmap(self.sites)
        for i in range(nsites):
            nbs = nbmap[i]
            if len(nbs) < 1:
                return False
        if not self.adsorbate.legit:
            return False
        if not self.defect.legit:
            return False

        def_ads_coords = self.defect.addsorption_coords(extend=0.0)[self.defect_adsid]
        # try:
        #     def_ads_coords = self.defect.addsorption_coords(extend=0.0)[self.defect_adsid]
        # except IndexError:
        #     print(self.defect.legit)
        #     print(self.defect.nolegit, print(self.defect.codename))
        #     self.defect.to('bug.xyz', 'xyz')
        #     self.to('config_bug.xyz', 'xyz')
        #     sys.exit('bug')
        ads_ads_coords = self.adsorbate.ads_site.coords
        if np.linalg.norm(ads_ads_coords - def_ads_coords) > 2.2 and self.adsorbate.speciesname != 'co':
            return False
        return True

    def cal_sg(self, mopwd, bincmd, mopheadline):
        out_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, self.codename, mopheadline,
                                                          caltype='sg')
        return out_dict, job_codename, mopin, mopout

    def cal_rlx(self, mopwd, bincmd, mopheadline):
        out_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, self.codename, mopheadline,
                                                          caltype='rlx')
        config_rlx = Config.from_outfile(mopout, self.defect.codename, self.adsorbate, self.defect_adsid,
                                         self.angle_ref, self.angle_perpen, caltype='rlx')
        return out_dict, job_codename, mopin, mopout, config_rlx

    def cal_thermal(self, mopwd, bincmd, mopheadline):
        out_dict, job_codename, mopin, mopout = mopac_run(self.sites, mopwd, bincmd, self.codename, mopheadline,
                                                          caltype='thermal')
        return out_dict, job_codename, mopin, mopout

    def to_mopin(self, mopwd, mopheadline, caltype):
        job_codename = self.codename + '-' + caltype
        mopin = '{}/{}.mop'.format(mopwd, job_codename)
        mopinstring = siteobjs2file(self.sites, fn=mopin, extension='mop', codename=job_codename,
                                    mopheadline=mopheadline)
        return mopinstring, mopin, job_codename


def out2out_dict(fn, codename, caltype):
    with open(fn, 'r') as f:
        ls = f.readlines()
    fstring = ''.join(ls)
    d = dict()
    fin = False
    sites = []
    tot = hof = gnorm = 0.0
    if caltype == 'rlx' or caltype == 'sg':
        igeo = 0
        for i in range(len(ls)):
            l = ls[i]
            if 'FINAL HEAT OF FORMATION' in l:
                fin = True
            if fin and 'FINAL HEAT OF FORMATION' in l:
                hof = float(l.strip().split()[5]) * 0.04336  # in eV
            elif fin and 'TOTAL ENERGY' in l:
                tot = float(l.strip().split()[3])
            elif fin and 'GRADIENT NORM' in l:
                gnorm = float(l.strip().split()[3])
            elif fin and 'CARTESIAN COORDINATES' in l:
                igeo = i
        if igeo != 0:
            for l in ls[igeo + 2:]:
                items = l.strip().split()
                if len(items) == 5:
                    num, e, x, y, z = items
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    sites.append(Site(e, [x, y, z]))
                else:
                    break

    elif caltype == 'thermal':
        ithermal = 0
        icoords = 0
        tot = None
        for i in range(len(ls)):
            l = ls[i]
            if '          CARTESIAN COORDINATES' in l:
                icoords = i
            if 'GRADIENT NORM = ' in l:
                gnorm = float(l.strip().split()[3])
            if 'CALCULATED THERMODYNAMIC PROPERTIES' in l:
                ithermal = i
                break
        for i in range(icoords+4, len(ls)):
            l = ls[i]
            items = l.strip().split()
            if len(items) == 5:
                num, e, x, y, z = items
                x = float(x)
                y = float(y)
                z = float(z)
                sites.append(Site(e, [x, y, z]))
            else:
                break
        if ithermal != 0:
            thermal_298_tot = ls[ithermal + 10].strip().split()[1:]
            hof, h, cp, entropy = thermal_298_tot
            hof = float(hof) * 0.043364
            entropy = float(entropy) * 0.043364 * 1e-3
            d['entropy'] = entropy
    d['tot'] = tot
    d['hof'] = hof
    d['sites'] = [s.as_dict() for s in sites]
    d['codename'] = codename
    d['fstring'] = fstring
    d['gnorm'] = gnorm
    return d


def mop2siteobjs(mopfn):
    with open(mopfn, 'r') as f:
        ls = f.readlines()
    ls = ls[3:]
    sites = []
    for l in ls:
        items = l.strip().split()
        if len(items) == 7:
            e, x, d1, y, d2, z, d3 = items
            x = float(x)
            y = float(y)
            z = float(z)
            sites.append(Site(e, [x, y, z]))
    return sites


def xyz2siteobjs(xyzfn):
    with open(xyzfn, 'r') as f:
        ls = f.readlines()
    ls = ls[2:]
    sites = []
    for l in ls:
        items = l.strip().split()
        if len(items) == 4:
            e, x, y, z = items
            x = float(x)
            y = float(y)
            z = float(z)
            sites.append(Site(e, [x, y, z]))
    return sites


def siteobjs2sitedicts(sites):
    return [s.as_dict() for s in sites]


def sitedicts2siteobjs(sitedicts):
    return [Site.from_dict(d) for d in sitedicts]


def siteobjs2string(sites, extension, codename, mopheadline=None):
    s = ''
    if extension == 'json':
        lst = siteobjs2sitedicts(sites)
        s += json.dumps(lst)

    elif extension == 'xyz':
        s += '{}\n{}\n'.format(len(sites), codename)
        for site in sites:
            s += '{}\t{:f} {:f} {:f}\n'.format(site.ele, site.x, site.y, site.z)

    elif extension == 'mop':
        if mopheadline is None:
            sys.exit('mop headline missing when export defect to mop file')
        s += '{}\n{}\n \n'.format(mopheadline, codename)
        for site in sites:
            s += '{}\t{:f} +1 {:f} +1 {:f} +1\n'.format(site.ele, site.x, site.y, site.z)
    return s


def siteobjs2file(sites, fn, extension, codename, mopheadline=None):
    s = siteobjs2string(sites, extension=extension, codename=codename, mopheadline=mopheadline)
    with open(fn, 'w') as f:
        f.write(s)
    return s


import os


def mopac_run(sites, mopwd, bincmd, codename, mopheadline, caltype):
    job_codename = codename + '-' + caltype
    mopin = '{}/{}.mop'.format(mopwd, job_codename)
    mopout = '{}/{}.out'.format(mopwd, job_codename)
    siteobjs2file(sites, fn=mopin, extension='mop', codename=job_codename, mopheadline=mopheadline)
    cmd = bincmd + ' ' + job_codename
    if os.path.isfile(mopout) and os.path.getsize(mopout) > 0:
        out_dict = out2out_dict(mopout, job_codename, caltype)
    else:
        process = subprocess.Popen(cmd, shell=True, cwd=mopwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        out_dict = out2out_dict(mopout, job_codename, caltype)
    return out_dict, job_codename, mopin, mopout


covalent_radii = {'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.73, 'N': 0.71, 'O': 0.66,
                  'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
                  'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
                  'Mn': 1.50, 'Fe': 1.42, 'Co': 1.38, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20,
                  'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
                  'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
                  'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 'Cs': 2.44, 'Ba': 2.15,
                  'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96,
                  'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
                  'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
                  'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21,
                  'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69}

vdw_radii = {'H': 1.2, 'He': 1.4, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47,
             'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.1, 'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Ar': 1.88,
             'K': 2.75, 'Ca': 2.31, 'Sc': 2.11, 'Ti': None, 'V': None, 'Cr': None, 'Mn': None, 'Fe': None,
             'Co': None, 'Ni': 1.63, 'Cu': 1.4, 'Zn': 1.39, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.9,
             'Br': 1.85, 'Kr': 0.88, 'Rb': 3.03, 'Sr': 2.49, 'Y': None, 'Zr': None, 'Nb': None, 'Mo': None,
             'Tc': None, 'Ru': None, 'Rh': None, 'Pd': 1.63, 'Ag': 1.72, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17,
             'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 1.08, 'Cs': 3.43, 'Ba': 2.68, 'La': None, 'Ce': None,
             'Pr': None, 'Nd': None, 'Pm': None, 'Sm': None, 'Eu': None, 'Gd': None, 'Tb': None, 'Dy': None,
             'Ho': None, 'Er': None, 'Tm': None, 'Yb': None, 'Lu': None, 'Hf': None, 'Ta': None, 'W': None,
             'Re': None, 'Os': None, 'Ir': None, 'Pt': 1.75, 'Au': 1.66, 'Hg': 1.55, 'Tl': 1.96, 'Pb': 2.02,
             'Bi': 2.07, 'Po': 1.97, 'At': 1.27, 'Rn': 1.2, 'Fr': None, 'Ra': None, 'Ac': None, 'Th': None,
             'Pa': None, 'U': None, 'Np': None, 'Pu': None, 'Am': None, 'Cm': None, 'Bk': None, 'Cf': None,
             'Es': None, 'Fm': None, 'Md': None, 'No': None, 'Lr': None, 'Rf': None, 'Db': None, 'Sg': None,
             'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None,
             'Mc': None, 'Lv': None, 'Ts': None, 'Og': None, }

atomic_numbers = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                  'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
                  'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
                  'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41,
                  'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
                  'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
                  'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69,
                  'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
                  'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88,
                  'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
                  'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107,
                  'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
                  'Ts': 117, 'Og': 118, }
