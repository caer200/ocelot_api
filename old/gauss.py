import fileop
import datetime
import multiprocessing
from pymatgen.core.structure import Molecule


class WriteInp:

    def __init__(self, jobtype, rawfile,
                 addkeywds='scf=(xqc,fermi,noincfock,ndamp=35,conver=8,vshift=500,novaracc) int=(acc2e=16)',
                 level='wb97xd/def2svp', mem='10GB',
                 nproc=multiprocessing.cpu_count()):
        self.jobtype = jobtype
        self.rawfile = rawfile
        self.addkeywds = addkeywds

        self.level = level
        self.mem = mem
        self.nproc = str(nproc)

        if self.jobtype == 'dimers':
            self.write_dimers()
        elif self.jobtype == 'reorg_opt':
            self.write_reorg_opt()
        elif self.jobtype == 'reorg_sp':
            self.write_reorg_sp()
        elif self.jobtype == 'tdsinglet':
            self.write_tdsinglet()

    def write_tdsinglet(self):
        WriteInp.writecom(jobname='tdsinglet',
                            nproc=self.nproc,
                            oldchk='nopt.chk',
                            mem=self.mem,
                            opt='td(NStates=20)',
                            level=self.level,
                            keywords=' geom=checkpoint ' + self.addkeywds,
                            charge=0,
                            spin=1,
                            sites=[])

    def write_dimers(self):
        raw = fileop.read_obj(self.rawfile)
        dimers = raw
        for dimer in dimers:
            self.dimer2com(dimer, self.level, self.addkeywds, self.nproc, self.mem)

    @staticmethod
    def dimer2com(dimer, level, addkeywds, nproc, mem):
        ref_mol = dimer.ref_organicmol
        oth_mol = dimer.oth_organicmol
        dimer_label = 'dimer' + dimer.label
        WriteInp.writecom(jobname=dimer_label + '_A',
                            nproc=nproc,
                            oldchk='',
                            mem=mem,
                            opt='',
                            level=level,
                            keywords='nosymm pop=full iop(3/33=1,3/59=8) int(acc2e=12) ' + addkeywds,
                            charge=0,
                            spin=1,
                            sites=ref_mol.mol.sites)
        WriteInp.writecom(jobname=dimer_label + '_B',
                            nproc=nproc,
                            oldchk='',
                            mem=mem,
                            opt='',
                            level=level,
                            keywords='nosymm pop=full iop(3/33=1,3/59=8) int(acc2e=12) ' + addkeywds,
                            charge=0,
                            spin=1,
                            sites=oth_mol.mol.sites)
        WriteInp.writecom(jobname=dimer_label + '_A_B',
                            nproc=nproc,
                            oldchk='',
                            mem=mem,
                            opt='',
                            level=level,
                            keywords='nosymm pop=full iop(3/33=1,3/59=8) int(acc2e=12) ' + addkeywds,
                            charge=0,
                            spin=1,
                            sites=ref_mol.mol.sites + oth_mol.mol.sites)

    def write_reorg_opt(self):
        raw = Molecule.from_file(self.rawfile)
        WriteInp.writecom(jobname='nopt',
                            nproc=self.nproc,
                            oldchk='',
                            mem=self.mem,
                            opt='opt',
                            level=self.level,
                            keywords='' + self.addkeywds,
                            charge=0,
                            spin=1,
                            sites=raw.sites)
        WriteInp.writecom(jobname='copt',
                            nproc=self.nproc,
                            oldchk='',
                            mem=self.mem,
                            opt='opt',
                            level=self.level,
                            keywords='' + self.addkeywds,
                            charge=1,
                            spin=2,
                            sites=raw.sites)
        WriteInp.writecom(jobname='aopt',
                            nproc=self.nproc,
                            oldchk='',
                            mem=self.mem,
                            opt='opt',
                            level=self.level,
                            keywords='' + self.addkeywds,
                            charge=-1,
                            spin=2,
                            sites=raw.sites)

    def write_reorg_sp(self):
        WriteInp.writecom(jobname='ncsp',
                            nproc=self.nproc,
                            oldchk='nopt.chk',
                            mem=self.mem,
                            opt='',
                            level=self.level,
                            keywords='geom=checkpoint' + self.addkeywds,
                            charge=1,
                            spin=2,
                            sites=[])
        WriteInp.writecom(jobname='nasp',
                            nproc=self.nproc,
                            oldchk='nopt.chk',
                            mem=self.mem,
                            opt='',
                            level=self.level,
                            keywords='geom=checkpoint' + self.addkeywds,
                            charge=-1,
                            spin=2,
                            sites=[])
        WriteInp.writecom(jobname='cnsp',
                            nproc=self.nproc,
                            oldchk='copt.chk',
                            mem=self.mem,
                            opt='',
                            level=self.level,
                            keywords='geom=checkpoint' + self.addkeywds,
                            charge=0,
                            spin=1,
                            sites=[])
        WriteInp.writecom(jobname='ansp',
                            nproc=self.nproc,
                            oldchk='aopt.chk',
                            mem=self.mem,
                            opt='',
                            level=self.level,
                            keywords='geom=checkpoint' + self.addkeywds,
                            charge=0,
                            spin=1,
                            sites=[])

    @staticmethod
    def writecom(jobname, nproc, oldchk, mem, level, opt, keywords, charge, spin, sites):
        """

        :param jobname:
        :param nproc:
        :param oldchk:
        :param mem:
        :param level:
        :param opt:
        :param keywords:
        :param charge:
        :param spin:
        :param sites:
        :return:
        """
        filename = jobname + '.com'
        with open(filename, 'w') as f:
            f.write('%rwf=' + jobname + '.rwf\r\n')
            f.write('%nosave\r\n')
            if oldchk != '':
                f.write('%oldchk=' + oldchk + '\r\n')
            f.write('%chk=' + jobname + '.chk\r\n')
            f.write('%nprocshared=' + str(nproc) + '\r\n')
            f.write('%mem=' + str(mem) + '\r\n')
            f.write('#p ' + opt + ' ' + level + ' ' + keywords + ' \r\n')
            f.write('\r\n')
            f.write('Input Generated at ' + str(datetime.datetime.now()) + ' for ' + jobname + '\r\n')
            f.write('\r\n')
            f.write(str(charge) + ' ' + str(spin) + ' \r\n')
            if len(sites) != 0:
                for s in sites:
                    x = str(round(s.x, 5))
                    y = str(round(s.y, 5))
                    z = str(round(s.z, 5))
                    f.write(s.species_string + '\t\t' + x + '  ' + y + '  ' + z + '  \r\n')
            else:
                f.write('\r\n')
            f.write('\r\n')
