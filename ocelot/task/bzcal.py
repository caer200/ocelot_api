from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import BSVasprun
from scipy import interpolate

from ocelot.routines import fileop
from ocelot.routines import mathop
from ocelot.routines.geometry import frac2cart
from ocelot.routines.geometry import norm

plt.switch_backend('agg')


class DispersionRelationLine:

    def __init__(self, scheme: str, reci_mat: np.ndarray, kcoords: np.ndarray, eigens_up: np.ndarray,
                 eigens_down: np.ndarray,
                 ivb: int, icb: int, linekdata=None, efermi=None):
        self.scheme = scheme
        self.reci_mat = reci_mat
        self.kcoords = kcoords
        self.ivb = ivb
        self.icb = icb
        self.eigens_up_occu = eigens_up  # eigens[ikpt][iband] = eigen, occu
        self.eigens_down_occu = eigens_down  # eigens[ikpt][iband] = eigen, occu, maybe None
        self.eigens_up = self.strip_occu_in_eigens(self.eigens_up_occu)
        if self.eigens_down_occu is None:
            self.eigens_down = None
        else:
            self.eigens_down = self.strip_occu_in_eigens(self.eigens_down_occu)

        self.linekdata = linekdata
        self.efermi = efermi

    @classmethod
    def from_vasprun_file(cls, vasprunfile):
        vdata = cls.read_vasprun(vasprunfile)
        return cls(
            'line',
            vdata['reci_mat'],
            vdata['kcoords'],
            vdata['eigens_up'],
            vdata['eigens_down'],
            ivb=int(int(vdata['nelect']) / 2 - 1),
            icb=int(int(vdata['nelect']) / 2),
            linekdata=None,
            efermi=vdata['efermi'],
        )



    def plotlinebs(self, fname='bands.eps', fitdensity=50, nexbands=2, zerofermi=True, channel='up',
                   yrange=(-2.0, 2.0)):
        f, ax = plt.subplots()
        fit_data = self.fit_linedata(fitdensity, nexbands, channel)
        if zerofermi:
            ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$'
            shift = self.efermi
        else:
            ylabel = r'$\mathrm{E\ (eV)}$'
            shift = 0.0
        for iband in range(self.ivb - nexbands, self.icb + nexbands + 1):
            ks_scatter = fit_data[iband]['ks']
            es_scatter = fit_data[iband]['es']
            ks_fit = fit_data[iband]['fitks']
            es_fit = fit_data[iband]['fites']
            if zerofermi:
                es_scatter -= shift
                es_fit -= shift
            ax.plot(ks_fit, es_fit, 'b-', linewidth=0.5)
            ax.scatter(ks_scatter, es_scatter, c='r', s=0.2)

        ax.set_xlabel(r'$\mathrm{k\ (\AA^{-1}})}$')
        ax.set_ylabel(ylabel)
        ax.set_ylim(yrange)
        ax.set_xlim((self.linekdata['xticks'][0], self.linekdata['xticks'][-1]))
        ax.set_xticks(self.linekdata['xticks'])
        ax.set_xticklabels(self.linekdata['xtickslabel'])
        for i in self.linekdata['xticks']:
            ax.axvline(x=i, color='black')
        plt.tight_layout()
        plt.savefig(fname, dpi=600)

    @staticmethod
    def strip_occu_in_eigens(eigens_with_occu: np.ndarray):
        size = np.shape(eigens_with_occu)[:2]
        eigens = np.zeros(size)
        ri, rj = size
        for i in range(ri):
            for j in range(rj):
                eigens[i, j] = eigens_with_occu[i, j, 0]
        return eigens

    def fit_linedata(self, fitdensity=50, nexbands=3, channel='up'):
        if channel == 'up':
            eigens = self.eigens_up
        else:
            if self.eigens_down is None:
                raise TypeError('down spin channel not available!')
            else:
                eigens = self.eigens_down
        fitted_data = {}
        for iband in range(self.ivb - nexbands, self.icb + nexbands + 1):
            fitted_data[iband] = {}
            fit_ks = []
            fit_es = []
            ks = []
            es = []
            lsegs = self.parse2linsegs(self.linekdata, iband, eigens, self.kcoords, self.reci_mat)
            for lseg in lsegs:
                seg_ks, seg_es, seg_rks, seg_res = lseg.fit(fitdensity)
                ks += seg_rks
                es += seg_res
                fit_ks += seg_ks
                fit_es += seg_es
            fitted_data[iband]['fitks'] = fit_ks
            fitted_data[iband]['fites'] = fit_es
            fitted_data[iband]['ks'] = ks
            fitted_data[iband]['es'] = es
            for k in fitted_data.keys():
                for kk in ['fitks', 'fites', 'ks', 'es']:
                    fitted_data[k][kk] = np.array(fitted_data[k][kk])
        return fitted_data

    @staticmethod
    def parse2linsegs(kdata, iband, eigens, kcoords, reci_mat):
        kline_density = kdata['kline_density']
        iseg = 0
        linsegs = []
        for seg_label in kdata['segs_label']:
            begin_label, end_label = seg_label.split('-')
            seg_kcoords = kcoords[iseg * kline_density: (iseg + 1) * kline_density]
            seg_eigens = np.array([e[iband] for e in eigens[iseg * kline_density: (iseg + 1) * kline_density]])
            iband = iband
            reci_mat = reci_mat
            pltstart = kdata['xticks'][iseg]
            lseg = LineSegment(begin_label, end_label, seg_kcoords, seg_eigens, iband, reci_mat, pltstart)
            linsegs.append(lseg)
            iseg += 1
        return linsegs

    @staticmethod
    def read_kpoints_line(kpointsfile, reci_mat):
        """
        kpoints must be a file generated by aflow

        return:

        d['kline_density'] = kline_density
        d['segs_label'] = segs_label
        d['xticks'] = kpath_kcoord
        d['xtickslabel'] = kpath_labels
        d['x'] = kpt_cumulative
        """
        with open(kpointsfile, 'r') as f_KPOINTS:
            lines = list(fileop.nonblank_lines(f_KPOINTS))

        kline_density = int(lines[1].split()[0])
        aflowstringkpath = [i.split('-') for i in lines[0].split()[2:]]
        segs_label = []
        for path in aflowstringkpath:
            if len(path) == 2:
                segs_label.append(path[0] + '-' + path[1])
            else:
                for i in range(len(path) - 1):
                    segs_label.append(path[i] + '-' + path[i + 1])

        kpath_kcoord = [0.0]  # these are xticks
        kroute = 0.0
        for i in range(0, len(lines[4:]), 2):
            words_s = lines[4:][i].split()
            words_e = lines[4:][i + 1].split()
            k_a_s = float(words_s[0])
            k_a_e = float(words_e[0])
            k_b_s = float(words_s[1])
            k_b_e = float(words_e[1])
            k_c_s = float(words_s[2])
            k_c_e = float(words_e[2])
            frac_v = np.array([k_a_e - k_a_s, k_b_e - k_b_s, k_c_e - k_c_s])
            length = np.linalg.norm(frac2cart(frac_v, reci_mat))
            kroute = kroute + length
            kpath_kcoord.append(kroute)

        kpath_labels = []  # these are xtickslabel
        for i in range(len(aflowstringkpath)):
            if i == 0:
                pre_seg_end = ''
            else:
                pre_seg_end = aflowstringkpath[i - 1][-1]
            for j in range(len(aflowstringkpath[i])):
                if j != 0 and j != len(aflowstringkpath[i]) - 1:
                    kpath_labels.append(aflowstringkpath[i][j])
                elif j == 0:
                    kpath_labels.append(pre_seg_end + '|' + aflowstringkpath[i][j])
        kpath_labels.append(aflowstringkpath[-1][-1])

        kpt_cumulative = np.zeros(((len(kpath_kcoord) - 1) * kline_density))  # these are x
        for i in range(len(kpath_kcoord) - 1):
            seg = np.linspace(kpath_kcoord[i], kpath_kcoord[i + 1], num=kline_density)
            for j in range(kline_density):
                kpt_cumulative[i * kline_density + j] = seg[j]
        d = OrderedDict()
        d['kline_density'] = kline_density
        d['segs_label'] = segs_label
        d['xticks'] = kpath_kcoord
        d['xtickslabel'] = kpath_labels
        d['x'] = kpt_cumulative
        return d

    @staticmethod
    def read_vasprun(vasprunfile):
        vasprun = BSVasprun(vasprunfile)
        vdata = vasprun.as_dict()
        # pprint(vdata['input']['crystal']['lattice']['matrix'])
        # pprint(vdata['input']['lattice_rec']['matrix'])
        real_latt = vdata['input']['crystal']['lattice']['matrix']
        real_latt = np.array(real_latt)
        reci_matrix = 2 * np.pi * np.linalg.inv(real_latt).T
        eigens_all = vasprun.eigenvalues
        if Spin(-1) not in eigens_all.keys():
            eigens_up = eigens_all[Spin(1)]
            eigens_down = None
        else:
            eigens_up = eigens_all[Spin(1)]
            eigens_down = eigens_all[Spin(-1)]  # eigens[ikpt][iband] = eigen, occu
        nelect = vdata['input']['parameters']['NELECT']
        efermi = vasprun.efermi
        kcoords = []
        for kpt in vdata['input']['kpoints']['actual_points']:
            kcoords.append(kpt['abc'])

        d = OrderedDict()
        d['reci_mat'] = reci_matrix
        d['nelect'] = int(nelect)
        d['efermi'] = efermi
        d['eigens_up'] = eigens_up
        d['eigens_down'] = eigens_down
        d['kcoords'] = np.array(kcoords)
        return d

    def get_line_ems(self, bandtype='vb', channel='up'):
        ems = []
        if bandtype == 'vb':
            iband = self.ivb
        elif bandtype == 'cb':
            iband = self.icb
        else:
            raise ValueError('bandtype should be cb or vb!')
        if channel == 'up':
            eigens = self.eigens_up
        elif channel == 'down':
            eigens = self.eigens_down
        else:
            raise ValueError('channel should be up or down!')
        if eigens is None:
            raise TypeError('eigens is None!')
        lsegs = self.parse2linsegs(self.linekdata, iband, eigens, self.kcoords, self.reci_mat)
        for lseg in lsegs:
            ems += lseg.get_ems(bandtype)
        for i in range(len(ems)):
            ipt, lseg, em = ems[i]
            ems[i] += [lseg.eigens[ipt] - self.efermi, lseg.dispersion]
        return ems


    @classmethod
    def from_klines_files(cls, vasprunfile, kpointsfile):
        vdata = cls.read_vasprun(vasprunfile)
        kdata = cls.read_kpoints_line(kpointsfile, vdata['reci_mat'])
        return cls(
            'line',
            vdata['reci_mat'],
            vdata['kcoords'],
            vdata['eigens_up'],
            vdata['eigens_down'],
            ivb=int(int(vdata['nelect']) / 2 - 1),
            icb=int(int(vdata['nelect']) / 2),
            linekdata=kdata,
            efermi=vdata['efermi'],
        )


class LineSegment:

    def __init__(self, begin_label: str, end_label: str, kcoords: np.ndarray, eigens: np.ndarray,
                 iband: int, reci_mat: np.ndarray, pltstart: float = 0.0):

        self.begin_label = begin_label
        self.end_label = end_label
        self.kcoords = kcoords
        self.linedensity = len(self.kcoords)
        self.eigens = eigens
        self.iband = iband
        self.reci_mat = reci_mat
        self.pltstart = pltstart  # the x coord in band structure plt for the first k point in this segment
        if len(kcoords) != len(eigens):
            raise ValueError('kcoords len != eigens len')

        self.kcumulative = np.linspace(self.pltstart, self.pltstart + self.length, self.linedensity)
        self.dispersion = np.amax(self.eigens) - np.amin(self.eigens)

    def fit(self, fitdensity=100):
        """
        fit band structure with spline at a predefined fit density
        using `scipy.interpolate.splev`

        :param int fitdensity: number of data points along one SEGMENT
        :return:
        """
        all_x = []
        all_y = []
        x_r = self.kcumulative
        y_r = self.eigens
        tck = interpolate.splrep(x_r, y_r, s=0.00)
        step = abs(x_r[0] - x_r[-1]) / fitdensity
        xs = [x * step + x_r[0] for x in range(fitdensity)]
        ys = [interpolate.splev(x * step + x_r[0], tck, der=0) for x in range(fitdensity)]
        ys_flat = list(np.array(ys).transpose())
        all_x += xs
        all_y += ys_flat
        return all_x, all_y, x_r.tolist(), y_r.tolist()

    def get_ems(self, bandtype='cb'):
        gmax, igmax = self.max_wi
        gmin, igmin = self.min_wi
        ems = []
        if bandtype == 'cb':
            poi = self.get_points_of_interest(self.eigens, gmin, bandtype, 0.025)
            for ipt in poi:
                em = self.get_em(self.kcumulative, self.eigens, ipt, 1)
                if em > 0.0:
                    ems.append([ipt, self, em])
        elif bandtype == 'vb':
            poi = self.get_points_of_interest(self.eigens, gmax, bandtype, 0.025)
            for ipt in poi:
                em = self.get_em(self.kcumulative, self.eigens, ipt, 1)
                if em < 0.0:
                    ems.append([ipt, self, em])
        return ems


    def __repr__(self):
        return "{}-{}:{}-{}".format(self.begin_label, self.end_label, self.kcoords[0], self.kcoords[-1])

    @property
    def pltend(self):
        return self.pltstart + self.length

    @property
    def frac_direction(self):
        return norm(self.kcoords[-1] - self.kcoords[0])

    @property
    def length(self):
        """
        length in k space
        """
        v = self.kcoords[-1] - self.kcoords[0]
        return norm(frac2cart(v, self.reci_mat))

    @property
    def max_wi(self):
        return np.amax(self.eigens), np.argmax(self.eigens)

    @property
    def min_wi(self):
        return np.amin(self.eigens), np.argmin(self.eigens)

    @staticmethod
    def get_em(kcumultive, eigens, ipt: int, step=1):
        rb_x = mathop.ra_to_rb(kcumultive)  # A-1 to b-1
        ha_y = mathop.ev_to_ha(eigens)
        em = mathop.fd_reci_2ndder(ha_y, rb_x, ipt, step=step)
        return em

    @staticmethod
    def get_points_of_interest(eigens, global_ex, bandtype: str, fluctuation=0.050):
        poi = []
        for j in range(len(eigens)):
            if bandtype == 'cb':
                if eigens[j] < global_ex + fluctuation:
                    poi.append(j)
                    # if (j == 0 and eigens[j] < eigens[j + step]) or (j == len(kcoords) - 1 and
                    #                                                  eigens[j] < eigens[j - step]):
                    #     poi.append(j)
                    # elif 2 * step < j < len(kcoords) - 2 * step and eigens[j] < eigens[j + step] and eigens[j] < eigens[
                    #     j - step]:
                    #     poi.append(j)
            elif bandtype == 'vb':
                if eigens[j] > global_ex - fluctuation:
                    poi.append(j)
                    # if (j == 0 and eigens[j] > eigens[j + step]) or (j == len(kcoords) - 1 and
                    #                                                  eigens[j] > eigens[j - step]):
                    #     poi.append(j)
                    # elif 2 * step < j < len(kcoords) - 2 * step and eigens[j] > eigens[j + step] and eigens[j] > eigens[
                    #     j - step]:
                    #     poi.append(j)
        return poi
