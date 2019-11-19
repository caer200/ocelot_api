import warnings
import os
from ocelot.routines.fileop import movefile
import subprocess
from pymatgen.core.structure import Element, Molecule
import numpy as np

# http://www.esi.umontreal.ca/accelrys/life/insight2000.1/zindo/3_Implementation.html
Zindo_elements = [
    'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K',
    'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Se', 'Br', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd'
]

ZINDO_INP_TEMPLATE = """
 $TITLEI
 
 {title}
 
 $END
 
 $CONTRL
 
 SCFTYP    {scftyp}      
 RUNTYP    {runtyp}        
 ENTTYP    {enttyp}
 UNITS     {units}     
 IPRINT    {iprint}
 INTTYP    {inttyp}
 IAPX      {iapx}
 MULT      {mult}   
 ITMAX     {itmax}   
 SCFTOL    {scftol}
 INTFA(1) = 1.00000 1.26700 0.58500 1.00000 1.00000 1.00000
 ONAME    = {oname}
 
 $END
 
 $DATAIN
 
 {datain}
 
 $END
"""

optional = """
 DYNAL(1) = {nzero} {ns} {nsp} {nspd} {nspdf} {cisize} {nactorb}
 NEL       {nel}

"""


class ZindoJob:

    def __init__(self, jobname, pmgmol, isdimer=False, mol_A=None, mol_D=None):
        self.jobname = jobname
        self.pmgmol = pmgmol  # add deepcopy?
        self.isdimer = isdimer
        self.mol_A = mol_A
        self.mol_D = mol_D

    @property
    def legit(self):
        if self.isdimer:
            sites = self.mol_A.sites + self.mol_D.sites
        else:
            sites = self.pmgmol.sites
        for s in sites:
            element = Element(s.species_string)
            if element.symbol not in Zindo_elements:
                return False
        return True

    @staticmethod
    def get_valence_electrons(pmgmol):
        """
        pymatgen has a weird `valence` property for element, e.g. for carbon valence is 2
        this is just for neutral mol

        :param pmgmol:
        :return:
        """
        nvelect = 0
        for site in pmgmol.sites:
            element = Element(site.species_string)
            nvelect += abs(min(element.common_oxidation_states))
        return nvelect

    @staticmethod
    def inpstring(pmgmol, RUNTYP='ENERGY', SCFTYP='RHF', ITMAX=500, SCFTOL=0.0000010, CISIZE=0, ACTSPC=0,
                  ONAME='zindoout'):
        """
        string of input file, charge and mult will be inherited from pmgmol

        :param pmgmol: pymatgen mol object
        :param RUNTYP: can be 'ENERGY' or 'GEO' or 'RPA' or 'CI'
        :param SCFTYP: type of scf
        :param ITMAX: max # of iteration
        :param SCFTOL: scf tol
        :param CISIZE: CI basis size
        :param ACTSPC: # of active orbitals in CI
        :param ONAME: output name
        :return:
        """
        if pmgmol.charge != 0 or pmgmol.spin_multiplicity != 1:
            warnings.warn('W: zindo io can only handle neutral singlet right now!')
        # nelect, ns, nsp, nspd, nspdf = [0]*5
        datain = ''
        for site in pmgmol.sites:
            element = Element(site.species_string)
            # nelect += abs(min(element.common_oxidation_states))
            datain += '{} {} {} {}\n'.format(site.x, site.y, site.z, element.number)
            # if element.block == 's':
            #     ns += 2
            # elif element.block == 'p':
            #     nsp += 2
            # elif element.block == 'd':
            #     nspd += 1
            # elif element.block == 'f':
            #     nspdf += 1
            # else:
            #     warnings.warn('W: this element has a row number > 7 ??')
        s = ZINDO_INP_TEMPLATE.format(
            title='zindogen', scftyp=SCFTYP, runtyp=RUNTYP, enttyp='COORD', units='ANGS', iprint='1', inttyp='1',
            iapx='3', mult=1, itmax=ITMAX, scftol=SCFTOL,
            # nel=nelect, nzero=0, ns=ns, nsp=nsp, nspd=nspd, nspdf=nspdf, cisize=CISIZE, nactorb=ACTSPC,
            oname=ONAME, datain=datain,
        )
        return s
    
    def write_single(self, fnprefix, RUNTYP='ENERGY', SCFTYP='RHF', ITMAX=500, SCFTOL=0.0000010, CISIZE=0, ACTSPC=0,
                     ONAME='zindoout'):
        s = self.inpstring(self.pmgmol, RUNTYP, SCFTYP, ITMAX, SCFTOL, CISIZE, ACTSPC, ONAME)
        with open(fnprefix + '.inp', 'w') as f:
            f.write(s)
        return [fnprefix + '.inp']

    def write_dimer(self, fnprefix, RUNTYP='ENERGY', SCFTYP='RHF', ITMAX=500, SCFTOL=0.0000010, CISIZE=0, ACTSPC=0,
                     ONAME='zindoout'):
        if not self.isdimer:
            return
        s_a = self.inpstring(self.mol_A, RUNTYP, SCFTYP, ITMAX, SCFTOL, CISIZE, ACTSPC, ONAME)
        s_d = self.inpstring(self.mol_D, RUNTYP, SCFTYP, ITMAX, SCFTOL, CISIZE, ACTSPC, ONAME)
        s_dimer = self.inpstring(self.pmgmol, RUNTYP, SCFTYP, ITMAX, SCFTOL, CISIZE, ACTSPC, ONAME)
        with open(fnprefix + '_dimer.inp', 'w') as f:
            f.write(s_dimer)
        with open(fnprefix + '_A.inp', 'w') as f:
            f.write(s_a)
        with open(fnprefix + '_D.inp', 'w') as f:
            f.write(s_d)
        return [
            fnprefix + '_dimer.inp',
            fnprefix + '_A.inp',
            fnprefix + '_D.inp',
        ]

    @classmethod
    def dimerjob_from_two_molecules(cls, pmgmol1, pmgmol2, jobname='dimer'):
        dimerpmgmol = Molecule.from_sites(pmgmol1.sites + pmgmol2.sites)
        return cls(jobname, dimerpmgmol, isdimer=True, mol_A=pmgmol1, mol_D=pmgmol2)

    def parse_tmo(self, tmofn='tmo.dat'):
        if not self.isdimer:
            warnings.warn('W: cannot parse tmo as this is not a dimer run!')
            return None

        nmo_a = self.get_valence_electrons(self.mol_A)
        nmo_d = self.get_valence_electrons(self.mol_D)
        data = np.empty((nmo_a, nmo_d), dtype=dict)
        with open(tmofn, 'r') as f:
            ls = f.readlines()

        for l in ls[3:]:
            imoa, imod, ti, ticm, emoa, emod = [float(number) for number in l.split()]
            imoa = int(imoa)
            imod = int(imod)
            d = dict(ti=ti, emoa=emoa, emod=emod, imoa=imoa, imod=imod)
            data[imoa-1][imod-1] = d
        return data

    @staticmethod
    def dimer_run(jobname, wdir, zindobin, zindoctbin, zindolib, pmgmola, pmgmolb):
        zj = ZindoJob.dimerjob_from_two_molecules(pmgmola, pmgmolb, jobname=jobname)
        whereami = os.getcwd()
        os.chdir(wdir)

        dimer_inps = zj.write_dimer(jobname)
        dimerinp, ainp, dinp = dimer_inps
        aout, dout, dimerout = ['A.out', 'D.out', 'dimer.out']
        env = dict(os.environ)
        try:
            ldlibpath = env['LD_LIBRARY_PATH']
        except KeyError:
            ldlibpath = ''
        env['LD_LIBRARY_PATH'] =  '{}:'.format(zindolib)+ ldlibpath
        p = subprocess.Popen("{} < {} > {}".format(zindobin, ainp, aout), shell=True, env=env)
        p.wait()
        movefile('mo.out', 'mo_A.out')
        p = subprocess.Popen("{} < {} > {}".format(zindobin, dinp, dout), shell=True, env=env)
        p.wait()
        movefile('mo.out', 'mo_D.out')
        p = subprocess.Popen("{} < {} > {}".format(zindoctbin, dimerinp, dimerout), shell=True, env=env)
        p.wait()
        data = zj.parse_tmo('tmo.dat')
        os.chdir(whereami)
        return data




