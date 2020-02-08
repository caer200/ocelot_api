import subprocess

from ocelot.routines.mopac import MopacInput, MopacOutput

"""
identify charged site with mopac routines

http://openmopac.net/Manual/charges.html
http://openmopac.net/Manual/Lewis_structures.html
"""


class IdChg:

    def __init__(self, pmgmol, jobname, mopaccmd):
        self.pmgmol = pmgmol
        self.jobname = jobname
        self.sites = pmgmol.sites
        self.mopaccmd = mopaccmd

        self.inputname = self.jobname + '.mop'
        self.outputname = self.jobname + '.out'

    def write_inp(self, wdir):
        min = MopacInput(self.sites, mopheader='CHARGES', comment_line='')
        min.write_mopinput(wdir + '/' + self.inputname)

    @staticmethod
    def parse_out(outfn):
        with open(outfn, 'r') as f:
            fstring = f.read()
        mout = MopacOutput(fstring=fstring, caltype='CHARGES')
        data = mout.parse_charges()
        return data

    def run(self, wdir):
        process = subprocess.Popen(self.mopaccmd, shell=True, cwd=wdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        d = self.parse_out(wdir + '/' + self.outputname)
        return d
