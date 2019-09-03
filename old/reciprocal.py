import warnings
import copy
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import fileop
import mathop
from scipy import interpolate
plt.switch_backend('agg')


class ReciprocalAnalyzer:
    """
    analysis for band structure

    firstly, parse eignvalues from OUTCAR KPOINTS EIGENVAL into a dict called band_data
    """

    def __init__(self, jtype, postshift=0.0, fitdensity=100, yrange=(-2, 2), usefitdata=True, zerofermi=True):
        """
        :param jtype:  a string used to name output fig
        :param postshift:  a float to shift the eigen spectrum
        :param fitdensity:  the # of points in a seg you want to fit
        :param yrange:  output fig yrange
        :param usefitdata:  do we use fit data to calculate the curvature
        :param zerofermi: do we zero fermi level
        """
        shouldistart = 0
        if os.path.isfile('OUTCAR') and os.path.isfile('KPOINTS') and os.path.isfile('EIGENVAL'):
            if fileop.check_outcar():
                shouldistart = 1

        if shouldistart == 0:
            warnings.warn('ReciprocalAnalyzer cannot find vasp outputs at ' + os.getcwd())
            sys.exit(1)

        self.zerofermi = zerofermi
        self.usefitdata = usefitdata
        self.jtype = jtype
        self.postshift = postshift
        self.fitdensity = fitdensity
        self.yrange = yrange

        self.r_d = self.eigen_parser()
        self.f_d = self.bandsfit(self.r_d, 100)
        self.cbems = self.get_ems('cb', 0.025, self.f_d, 3)
        self.vbems = self.get_ems('vb', 0.025, self.f_d, 3)
        self.e_em_entry = sorted(self.cbems, key=lambda x: abs(x[0]))[0]
        self.h_em_entry = sorted(self.vbems, key=lambda x: abs(x[0]))[0]
        self.e_em = round(self.e_em_entry[0], 2)
        self.h_em = round(self.h_em_entry[0], 2)
        self.bands_plt(self.f_d, self.r_d, self.yrange, self.jtype)
        self.gap = min(self.f_d['y'][self.f_d['cb_idx']]) - max(self.f_d['y'][self.f_d['vb_idx']])

    @staticmethod
    def bands_plt(f_d, r_d, yrange, jtype):

        savefigname = 'vaspbands_' + jtype + '.png'

        f, ax = plt.subplots()

        for i in range(len(f_d['y'])):
            x_fit = f_d['x']
            y_fit = f_d['y'][i]
            ax.plot(x_fit, y_fit, 'b-', linewidth=0.5, alpha=0.9)

        for i in range(len(r_d['y'])):
            ax.scatter(r_d['x'], r_d['y'][i], c='r', s=0.2)

        ax.set_xlabel(r'$\mathrm{k\ (2\pi \cdot \AA^{-1}})}$')
        ylabel = r'$\mathrm{E\ -\ E_f\ (eV)}$'
        ax.set_ylabel(ylabel)
        ax.set_xlim(r_d['x'][0], r_d['x'][-1])
        ax.set_ylim(yrange)
        ax.set_xticks(r_d['xticks'])
        ax.set_xticklabels(r_d['xtickslabel'])
        for i in r_d['xticks']:
            ax.axvline(x=i, color='black')
        plt.tight_layout()
        plt.savefig(savefigname, dpi=600)

    @staticmethod
    def get_ems(band_type, fluctuation, data, step):
        """
        get ems for a certian band type
        :param band_type: 
        :param fluctuation: 
        :param data: 
        :param step: the step used to define the boundary of parabolla
        :return: 
        """
        kline_density = data['kline_density']
        segs_label = data['segs_label']
        x = data['x']
        if band_type == 'cb':
            y = data['y'][data['cb_idx']]
            g_ex = min(y)
            ems = []
            for i in range(len(segs_label)):
                seg_x = x[i * kline_density: (i + 1) * kline_density]
                seg_y = y[i * kline_density: (i + 1) * kline_density]
                poi = []  # point of interest
                for j in range(kline_density):
                    if (j == 0 and seg_y[j] < seg_y[j + step]) or (j == kline_density - 1 and
                                                                   seg_y[j] < seg_y[j - step]):
                        if seg_y[j] < g_ex + fluctuation:
                            poi.append(j)
                    elif 2 * step < j < kline_density - 2 * step:
                        if seg_y[j] < seg_y[j + step] and seg_y[j] < seg_y[j - step] and seg_y[j] < g_ex + fluctuation:
                            poi.append(j)
                for idx in poi:
                    rb_x = mathop.ra_to_rb(seg_x)
                    ha_y = mathop.ev_to_ha(seg_y)
                    em = mathop.fd_reci_2ndder(ha_y, rb_x, idx, step=step)
                    if em > 0.0:
                        ems.append((em, idx / (kline_density - 1), segs_label[i], band_type))
            return ems

        elif band_type == 'vb':
            y = data['y'][data['vb_idx']]
            g_ex = max(y)
            ems = []
            for i in range(len(segs_label)):
                seg_x = x[i * kline_density: (i + 1) * kline_density]
                seg_y = y[i * kline_density: (i + 1) * kline_density]
                poi = []  # point of interest
                for j in range(kline_density):
                    if (j == 0 and seg_y[j + step] < seg_y[j]) or (j == kline_density - 1 and
                                                                   seg_y[j - step] < seg_y[j]):
                        if seg_y[j] > g_ex - fluctuation:
                            poi.append(j)
                    elif 2 * step < j < kline_density - 2 * step:
                        if seg_y[j] > seg_y[j + step] and seg_y[j] > seg_y[j - step] and seg_y[j] > g_ex - fluctuation:
                            poi.append(j)
                for idx in poi:
                    rb_x = mathop.ra_to_rb(seg_x)
                    ha_y = mathop.ev_to_ha(seg_y)
                    em = mathop.fd_reci_2ndder(ha_y, rb_x, idx, step=step)
                    if em < 0.0:
                        ems.append((em, idx / (kline_density - 1), segs_label[i], band_type))
            return ems

    @staticmethod
    def bandsfit(data, fitdensity, nbandsud=5):

        d = copy.deepcopy(data)

        segs_label = data['segs_label']
        x_raw = data['x']
        kline_density = int(data['kline_density'])
        eigens = data['y']
        nelect = data['nelect']

        eigenstofit = []  # proxy bands 5+5 x nkpts
        for i in range(int((nelect / 2) - nbandsud), int((nelect / 2) + nbandsud)):
            eigenstofit.append(eigens[i])
        d['vb_idx'] = int(len(eigenstofit) / 2) - 1
        d['cb_idx'] = int(len(eigenstofit) / 2)

        fit_eigens = np.zeros((len(eigenstofit), fitdensity * len(segs_label)))

        all_x = []
        for iband in range(len(eigenstofit)):
            y_raw = eigenstofit[iband]
            all_x = []
            all_y = []
            for seg_index in range(len(segs_label)):
                x_r = x_raw[seg_index * kline_density: (seg_index + 1) * kline_density]
                y_r = y_raw[seg_index * kline_density: (seg_index + 1) * kline_density]
                tck = interpolate.splrep(x_r, y_r, s=0.00)
                step = abs(x_r[0] - x_r[-1]) / fitdensity
                xs = [x * step + x_r[0] for x in range(fitdensity)]
                ys = [interpolate.splev(x * step + x_r[0], tck, der=0) for x in range(fitdensity)]
                ys_flat = list(np.array(ys).transpose())
                all_x += xs
                all_y += ys_flat
            fit_eigens[iband] = all_y
        d['kline_density'] = fitdensity
        d['y'] = fit_eigens
        d['x'] = all_x
        return d

    @staticmethod
    def eigen_parser(zerofermi=True, postshift=0.0):
        """
        parse vasp output for eigenvalue problem
        :param zerofermi:
        :param postshift:
        :return: a dictionary
        """
        d = dict()

        r_mat = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]
        efermi = 0.0
        nelect = 0
        with open('OUTCAR', 'r') as f_OUTCAR:
            lines = f_OUTCAR.readlines()
        for i, line in enumerate(lines):
            if line.startswith('      direct lattice vectors'):
                v1 = [float(num) for num in lines[i+1].split()[3:]]
                v2 = [float(num) for num in lines[i+2].split()[3:]]
                v3 = [float(num) for num in lines[i+3].split()[3:]]
                r_mat = [v1, v2, v3]
            if line.startswith(' E-fermi'):
                efermi = float(lines[i].split()[2])
            if line.startswith('   NELECT'):
                nelect = int(round(float(lines[i].split()[2])))
        d['rmat'] = r_mat
        d['efermi'] = efermi
        d['nelect'] = nelect
        d['vb_idx'] = int(nelect / 2) - 1
        d['cb_idx'] = int(nelect / 2)

        with open('KPOINTS', 'r') as f_KPOINTS:
            lines = list(fileop.nonblank_lines(f_KPOINTS))

        kline_density = int(lines[1].split()[0])
        aflowstringkpath = [i.split('-') for i in lines[0].split()[2:]]
        segs_label = []
        for path in aflowstringkpath:
            if len(path) == 2:
                segs_label.append(path[0] + '-' + path[1])
            else:
                for i in range(len(path)-1):
                    segs_label.append(path[i] + '-' + path[i + 1])
        d['kline_density'] = kline_density
        d['segs_label'] = segs_label

        kpath_kcoord = [0.0]  # these are xticks
        kroute = 0.0
        for i in range(0, len(lines[4:]), 2):
            words_s = lines[4:][i].split()
            words_e = lines[4:][i+1].split()
            k_a_s = float(words_s[0])
            k_a_e = float(words_e[0])
            k_b_s = float(words_s[1])
            k_b_e = float(words_e[1])
            k_c_s = float(words_s[2])
            k_c_e = float(words_e[2])
            frac_v = np.array([k_a_e-k_a_s, k_b_e-k_b_s, k_c_e-k_c_s])
            length = np.linalg.norm(mathop.coords_converter2cart(frac_v, r_mat))
            kroute = kroute + length
            kpath_kcoord.append(kroute)
        d['xticks'] = kpath_kcoord

        kpath_labels = []  # these are xtickslabel
        for i in range(len(aflowstringkpath)):
            if i == 0:
                pre_seg_end = ''
            else:
                pre_seg_end = aflowstringkpath[i-1][-1]
            for j in range(len(aflowstringkpath[i])):
                if j != 0 and j != len(aflowstringkpath[i])-1:
                    kpath_labels.append(aflowstringkpath[i][j])
                elif j == 0:
                    kpath_labels.append(pre_seg_end+'|'+aflowstringkpath[i][j])
        kpath_labels.append(aflowstringkpath[-1][-1])
        d['xtickslabel'] = kpath_labels

        kpt_cumulative = np.zeros(((len(kpath_kcoord) - 1) * kline_density))  # these are x
        for i in range(len(kpath_kcoord) - 1):
            seg = np.linspace(kpath_kcoord[i], kpath_kcoord[i + 1], num=kline_density)
            for j in range(kline_density):
                kpt_cumulative[i * kline_density + j] = seg[j]
        d['x'] = kpt_cumulative

        with open('EIGENVAL', 'r') as f_EIGENVAL:
            lines = list(fileop.nonblank_lines(f_EIGENVAL))
        nkpts = int(lines[5].split()[1])
        nbands = int(lines[5].split()[2])
        eigens = np.empty([nbands, nkpts], dtype='float')  # these are y
        # kpts = []
        eigen_lines = lines[6:]
        for i in range(nkpts):
            block = eigen_lines[i*(nbands+1):(i+1)*(nbands+1)]
            # kpt_fcoord = [float(coord) for coord in block[0].split()[0:3]]
            # kpt_v = ops.geomath.coords_converter2cart(kpt_fcoord, r_mat)
            # kpts.append(kpt_v)
            for j in range(nbands):
                if zerofermi:
                    eigens[j][i] = float(block[j+1].split()[1]) - efermi - postshift
                else:
                    eigens[j][i] = float(block[j+1].split()[1]) - postshift
        d['y'] = eigens
        d['nkpts'] = nkpts
        d['nbands'] = nbands

        return d
