from ccdc.io import EntryReader
import sys
'''
these are the most general rules that descide whether a certain structure may be imported into fomdb from csd:
1. at least 2 fused rings where for each ring, there're at least 2 unsaturated bonds
2. there's no metal/open-shell atom in the system
3. no polymer
4. one kind of molecule

output:
1. write a csv, csd_id & csd_SMILES
2. write csd_id.cif
'''
csd_entry_reader = EntryReader('CSD')


def is_condition4(entry):
    if len(entry.molecule.components) == 1:
        return 1
    return 0


def is_condition3(entry):
    if not entry.molecule.components[0].is_polymeric:
        return 1
    return 0


def is_condition2(entry):
    atoms = entry.molecule.components[0].atoms
    if len([a for a in atoms if a.is_metal]) == 0:
        return 1
    return 0


def is_condition1(entry):
    frings = [r for r in entry.molecule.components[0].rings if r.is_fused and is_atleast2unsaturated(r)]
    if len(frings) > 1:
        return 1
    return 0


def is_atleast2unsaturated(ring):
    ubs = [b for b in ring.bonds if b.bond_type in (2, 3, 4, 5, 7, 9)]
    if len(ubs) > 1:
        return 1
    return 0


def roundfloat(f):
    if isinstance(f, float):
        return round(f, 2)
    return f


def main():
    column_idx = ['csd_id', 'csd_smiles']
    rescsv = open('csdcrawl.csv', 'a')
    rescsv.write(";;;" + ";;;".join(column_idx) + "\n")
    # main loop
    idxcounter = 0
    # for entry in [csd_entry_reader[i] for i in range(100)]:  # debug before launching
    for entry in csd_entry_reader:
        # add filters
        if not is_condition4(entry):
            continue
        if not is_condition3(entry):
            continue
        if not is_condition2(entry):
            continue
        if not is_condition1(entry):
            continue
        csdid = entry.identifier
        try:
            smiles = entry.molecule.smiles
        except RuntimeError:
            smiles = None
        rescsv.write(str(idxcounter) + ";;;" + ";;;".join([str(csdid), str(smiles)]) + "\n")
        cifstring = entry.to_string(format='cif')
        with open(str(csdid)+'.cif', 'w') as f:
            f.write(cifstring)
        idxcounter += 1
    rescsv.close()


def cleanup():
    sys.exit(1)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        cleanup()


