from ocelot.schema.msite import MSite
import numpy as np


class MopacInput:

    def __init__(self, msites, mopheader='CHARGES', comment_line=''):
        """
        does not support selective dynamics

        :param msites:
        :param mopheader:
        :param comment_line:
        """
        self.msites = msites
        self.mopheader = mopheader
        self.comment_line = comment_line

    @staticmethod
    def sites2mopinput(msites, header, comment):
        """
        return the string of mopac input

        :param sites:
        :return:
        """
        s = ""
        s += '{}\n{}\n \n'.format(header, comment)
        for msite in msites:
            s += '{}\t{:f} +1 {:f} +1 {:f} +1\n'.format(msite.element.name, msite.x, msite.y, msite.z)
        return s

    def write_mopinput(self, mopfn):
        s = self.sites2mopinput(self.msites, self.mopheader, self.comment_line)
        with open(mopfn, 'w') as f:
            f.write(s)

    @classmethod
    def from_mopinput(cls, mopfn):
        with open(mopfn, 'r') as f:
            ls = f.readlines()
        mopheader = ls[0].strip()
        comment_line = ls[1].strip()
        ls = ls[3:]
        sites = []
        for l in ls:
            items = l.strip().split()
            if len(items) == 7:
                e, x, d1, y, d2, z, d3 = items
                x = float(x)
                y = float(y)
                z = float(z)
                sites.append(MSite(e, [x, y, z]))
        return cls(sites, mopheader=mopheader, comment_line=comment_line)


class MopacOutput:

    def __init__(self, fstring, caltype):
        self.fstring = fstring
        self.caltype = caltype

    def parse_thermo(self):
        ls = self.fstring.split('\n')
        fin = False
        msites = []
        hof = gnorm = entropy = None
        ithermal = 0
        icoords = 0
        for i in range(len(ls)):
            l = ls[i]
            if '          CARTESIAN COORDINATES' in l:
                icoords = i
            if 'GRADIENT NORM = ' in l:
                gnorm = float(l.strip().split()[3])
            if 'CALCULATED THERMODYNAMIC PROPERTIES' in l:
                ithermal = i
            if '== MOPAC DONE ==' in l:
                fin = True
        for i in range(icoords + 4, len(ls)):
            l = ls[i]
            items = l.strip().split()
            if len(items) == 5:
                num, e, x, y, z = items
                x = float(x)
                y = float(y)
                z = float(z)
                msites.append(MSite(e, [x, y, z]))
            else:
                break
        if ithermal != 0:
            thermal_298_tot = ls[ithermal + 10].strip().split()[1:]
            hof, h, cp, entropy = thermal_298_tot
            hof = float(hof) * 0.043364
            entropy = float(entropy) * 0.043364 * 1e-3
        if fin:
            d = dict()
            d['caltype'] = 'thermo'
            d['gnorm'] = gnorm
            d['heat_of_formation_298'] = hof
            d['entropy'] = entropy
            d['msites'] = [s.as_dict() for s in msites]
            d['fstring'] = self.fstring
            return d
        return None

    def parse_rlx_or_sp(self):
        ls = self.fstring.split('\n')
        fin = False
        msites = []
        igeo = 0
        tot = hof = gnorm = None
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
                    msites.append(MSite(e, [x, y, z]))
                else:
                    break
        if fin:
            d = dict()
            d['caltype'] = 'thermo'
            d['gnorm'] = gnorm
            d['heat_of_formation'] = hof
            d['total_energy'] = tot
            d['msites'] = [s.as_dict() for s in msites]
            d['fstring'] = self.fstring
            return d
        return None

    def parse_charges(self):
        ls = self.fstring.split('\n')
        fin = False
        msites = []
        icharge = 0
        igeo = 0
        total_charge = 0
        for i in range(len(ls)):
            l = ls[i]
            if 'JOB ENDED NORMALLY' in l:
                fin = True
            if 'COMPUTED CHARGE ON SYSTEM:' in l:
                total_charge = int(l.strip().split()[-1])
            if '          CARTESIAN COORDINATES' in l:
                igeo = i
            if 'Ion Atom No.  Type    Charge' in l:
                icharge = i

        if igeo != 0:
            for l in ls[igeo + 4:]:
                items = l.strip().split()
                if len(items) == 5:
                    num, e, x, y, z = items
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    msites.append(MSite(e, [x, y, z]))
                else:
                    break
            charges = np.zeros(len(msites))
            for l in ls[icharge + 2:]:
                items = l.strip().split()
                if len(items) == 4:
                    iatom, atom_number, atom_type, atom_charge = items
                    iatom = int(iatom)
                    atom_charge = int(atom_charge)
                    charges[iatom] = atom_charge
                else:
                    break
            if fin:
                d = dict()
                d['total_charge'] = total_charge
                d['charges'] = charges
                return d
        return None


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
            sites.append(MSite(e, [x, y, z]))
    return sites
