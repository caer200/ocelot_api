from pymatgen.command_line.bader_caller import BaderAnalysis

"""
http://theory.cm.utexas.edu/henkelman/code/bader/
https://cms.mpi.univie.ac.at/wiki/index.php/LAECHG

One major issue with the charge density (CHGCAR) files from the VASP code is that they only contain the 
valance charge density. The Bader analysis assumes that charge density maxima are located at atomic centers 
(or at pseudoatoms). Aggressive pseudopotentials remove charge from atomic centers where it is both expensive 
to calculate and irrelevant for the important bonding properties of atoms.

Recently, the VASP developers have added a module (aedens) which allows for the core charge to be written out 
from PAW calculations. This module is included in vasp version 4.6.31 08Feb07 and later. 
By adding the LAECHG=.TRUE. to the INCAR file, the core charge is written to AECCAR0 and the valance charge to AECCAR2. 
These two charge density files can be summed using the chgsum.pl script:

  chgsum.pl AECCAR0 AECCAR2

The total charge will be written to CHGCAR_sum.

The bader analysis can then be done on this total charge density file:

  bader CHGCAR -ref CHGCAR_sum

One finally note is that you need a fine fft grid to accurately reproduce the correct total core charge. 
It is essential to do a few calculations, increasing NG(X,Y,Z)F until the total charge is correct. 
"""

ba = BaderAnalysis('CHGCAR', 'POTCAR')  # this is a good example for running code and capturing stdout
# print(ba.get_charge(1))

from ocelot.routines.pbc import PBCparser
from pymatgen.core.structure import Structure

struct = Structure.from_file('CONTCAR')
mols, unwrap_str_sorted, unwrap_pblock_list = PBCparser.unwrap(struct)


def parse_acf(acfdata_fn='./ACF.dat'):
    acfdata = []
    with open(acfdata_fn, 'r') as f:
        ls = f.readlines()
    for l in ls[2:-4]:
        iatom, x, y, z, charge, mindist, atomvol = l.strip().split()
        entry = {}
        entry['iatom'] = int(iatom) - 1
        entry['x'] = x
        entry['y'] = y
        entry['z'] = z
        entry['charge'] = charge
        entry['mindist'] = mindist
        entry['atomvol'] = atomvol
        for k in entry.keys():
            entry[k] = float(entry[k])
        acfdata.append(entry)
    return acfdata


# acfdata = parse_acf()
imols = [s.properties['imol'] for s in unwrap_str_sorted.sites]
isites = [s.properties['siteid'] for s in unwrap_str_sorted.sites]
istie2imol_dict = dict(zip(isites, imols))
for imol in range(len(mols)):
    molchg = 0
    for isite in range(len(isites)):
        if istie2imol_dict[isite] == imol:
            # molchg += acfdata[isite]['charge']
            molchg += ba.get_charge(isite)
    print(molchg)
