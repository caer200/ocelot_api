import os


class Reorganization:
    """ class to calculate reorganization energy for holes and electrons from Gaussian output files."""

    def __init__(self, nopt='nopt.log', copt='copt.log', aopt='aopt.log', ngeo_a='nasp.log', ngeo_c='ncsp.log',
                 ageo_n='ansp.log', cgeo_n='cnsp.log'):

        self.nopt = nopt
        self.copt = copt
        self.aopt = aopt

        self.ngeo_a = ngeo_a  # n geometry, charge spin as a
        self.ngeo_c = ngeo_c
        self.ageo_n = ageo_n
        self.cgeo_n = cgeo_n

        self.fns = (self.nopt, self.copt, self.aopt, self.ngeo_a, self.ngeo_c, self.ageo_n, self.cgeo_n)

        for fn in self.fns:
            if not (os.path.isfile(fn) and os.stat(fn).st_size > 0):
                raise FileNotFoundError('log files for reorg calculations not ready')

    def reorg_energy(self):

        aopt_energy = self.extract_energy(self.aopt)
        copt_energy = self.extract_energy(self.copt)
        nopt_energy = self.extract_energy(self.nopt)

        ngeo_a_energy = self.extract_energy(self.ngeo_a)
        ngeo_c_energy = self.extract_energy(self.ngeo_c)
        ageo_n_energy = self.extract_energy(self.ageo_n)
        cgeo_n_energy = self.extract_energy(self.cgeo_n)

        hole_reorg = (ngeo_c_energy - copt_energy + cgeo_n_energy - nopt_energy) * 27.2114
        elec_reorg = (ngeo_a_energy - aopt_energy + ageo_n_energy - nopt_energy) * 27.2114  # ev

        ip_adia = (ngeo_c_energy - nopt_energy) * 27.2114
        ip_nadia = (copt_energy - nopt_energy) * 27.2114

        ea_adia = (ngeo_a_energy - nopt_energy) * 27.2114
        ea_nadia = (aopt_energy - nopt_energy) * 27.2114

        return hole_reorg, elec_reorg, ip_adia, ip_nadia, ea_adia, ea_nadia

    @staticmethod
    def extract_energy(log):
        with open(log, 'r') as f:
            loglines = f.readlines()
        scflines = []
        mp2line = 0
        for line in loglines:
            if "SCF Done" in line:
                scflines.append(line.strip())
            elif "EUMP2" in line:
                mp2line = line.strip()
            elif "ERMP2" in line:
                mp2line = line.strip()
        if mp2line != 0:
            mp2line.replace('D', 'E')
            mp2line = mp2line.split()
            return float(mp2line[5])
        else:
            linelen = len(scflines)
            scfline = scflines[linelen - 1].strip().split()
            return float(scfline[4])
