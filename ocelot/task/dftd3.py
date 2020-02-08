import subprocess
import warnings
from ocelot.routines.fileop import createdir

"""
dftd3 from Grimme for pbc system

-func <functional name in TM style>
Choose one of the implemented functionals. No default. For a list of parameterized func-
tionals, refer to our web-site.9 HF can be invoked by -func hf.

- '-anal', '-grad', '-tz', options are not supported currently
- citations
    main code: S. Grimme, J. Antony, S. Ehrlich and H. Krieg, J. Chem. Phys. 132 (2010), 154104
    BJ-damping: S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem. 32 (2011), 1456-1465
    DFT-D2: S. Grimme, J. Comput. Chem., 27 (2006), 1787-1799
    DFT-D3M/DFT-D3M(BJ): D.G.A. Smith, L.A. Burns, K. Patkowski and C.D. Sherrill, J. Phys. Chem. Lett. 7 (2016) 2197-2203
"""
DAMPING_OPTIONS = ['bj', 'bjm', 'zero', 'zerom']
FUNCTIONAL_OPTIONS = [
    'HF', 'B1B95', 'B2GPPLYP', 'B3PW91', 'BHLYP', 'BMK', 'BOP', 'BPBE', 'CAMB3LYP', 'LCωPBE', 'MPW1B95', 'MPWB1K',
    'mPWLYP', 'OLYP', 'OPBE', 'oTPSS', 'PBE38', 'PBEsol', 'PTPSS', 'PWB6K', 'revSSB', 'SSB', 'TPSSh', 'HCTH120',
    'B2PLYP', 'B3LYP', 'B97D', 'BLYP', 'BP86', 'DSDBLYP', 'PBE0', 'PBE', 'PW6B95', 'PWPB95', 'revPBE0', 'revPBE38',
    'revPBE', 'TPSS0', 'TPSS',
]
FUNCTIONAL_OPTIONS = [func.lower() for func in FUNCTIONAL_OPTIONS]


class DFTD3:

    def __init__(self, jobname, structure, func='pbe', damping='bj', dftd2=False, cutoff=94.8683, cnthr=40):
        """
        init a dftd3 calculation

        :param jobname:
        :param structure:
        :param func:
        :param damping:
        :param dftd2: DFT-D2 version.10
        :param cutoff: a cutoff value for the dispersion interaction. The default value is 95 a.u.
        :param cnthr: a cutoff value for the calculation of the CN. The default value is 40 a.u. and should be kept fixed
        """
        self.jobname = jobname
        self.structure = structure
        if func not in FUNCTIONAL_OPTIONS:
            warnings.warn('W: the func you spec is {} which is not supported, im using pbe instead'.format(func))
            self.func = 'pbe'
        else:
            self.func = func
        if damping not in DAMPING_OPTIONS:
            warnings.warn('W: the damping you spec is {} which is not supported, im using bj instead'.format(damping))
            self.damping = 'bj'
        else:
            self.damping = damping
        self.dftd2 = dftd2
        self.cutoff = cutoff
        self.cnthr = cnthr

    @property
    def cmd_option_string(self):
        s = ' -func {} -{} -abc -pbc -cutoff {} -cnthr {}'.format(self.func, self.damping, self.cutoff, self.cnthr)
        if self.dftd2:
            s += ' -old'
        return s

    @staticmethod
    def parse_results(res):
        ls = res.split('\n')
        for l in ls:
            if 's6       :' in l:
                s6 = l.strip().split()[2]
                s6 = float(s6)
            elif 's8       :' in l:
                s8 = l.strip().split()[2]
                s8 = float(s8)
            elif 'a1       :' in l:
                a1 = l.strip().split()[2]
                a1 = float(a1)
            elif 'a2       :' in l:
                a2 = l.strip().split()[2]
                a2 = float(a2)
            elif 'Edisp /kcal,au,eV:' in l:
                edisp_ev = l.strip().split()[-1]
                edisp_ev = float(edisp_ev)
            elif 'E6    /kcal' in l:
                e6_kcal = l.strip().split()[-1]
                e6_ev = float(e6_kcal) * 0.04336
            elif 'E8    /kcal' in l:
                e8_kcal = l.strip().split()[-1]
                e8_ev = float(e8_kcal) * 0.04336
            elif 'E6(ABC)' in l:
                e6abc_kcal = l.strip().split()[-1]
                e6abc_ev = float(e6abc_kcal) * 0.04336
        try:
            d = dict(s6=s6, s8=s8, a1=a1, a2=a2, edisp_ev=edisp_ev, e6_ev=e6_ev, e6abc_ev=e6abc_ev, e8_ev=e8_ev)
        except NameError:
            return None
        return d

    def run(self, d3cmd, wdir):
        createdir(wdir)
        poscarfn = wdir + 'POSCAR_dftd3'
        self.structure.to('poscar', poscarfn)
        result = subprocess.check_output('{} {} {}'.format(d3cmd, poscarfn, self.cmd_option_string), shell=True)
        return result.decode('utf-8')
