from api.routines.pbc import CIFparser
from api.schema.config import Config
import glob
from api.task.packing import PackingIdentifier
import json

def get_tips():
    with open ("../../cifs/128.cif", "r") as myfile:
        data=myfile.read()
    cp = CIFparser.from_cifstring(data)
    cp.write_clean_cifs(prefix='tips')
    c = Config.from_file('{}-0.cif'.format('tips'))
    m = c.omols[0]
    tips_sc = sorted(m.sidechains, key=lambda x: len(x), reverse=True)[0]
    return tips_sc


def parser(cif, tips_sc):
    name = cif.split('/')[-1][:-4]
    with open (cif, "r") as myfile:
        data=myfile.read()
    cp = CIFparser.from_cifstring(data)
    cp.write_clean_cifs(prefix=name)
    c = Config.from_file('{}-0.cif'.format(name))
    c.unwrap_structure.to('cif', 'uconfig-{}.cif'.format(name))
    bc, bcpstructure, terminated_backbones= c.get_bone_config()
    bc.unwrap_structure.to('cif', 'bconfig-{}.cif'.format(name))
    pid = PackingIdentifier(bc)
    packingdata = pid.identify_heuristaic()
    m = [m for m in c.omols if not m.is_solvent][0]
    bone = m.backbone
    tips_sc_u = tips_sc.umbrella
    bd = {
        'nrings': len(bone.backbone_rings),
        'linearity': bone.lfit_linearity,
        'lp': bone.lp,
        'lq': bone.lq,
        'canonical_smiles': bone.canonical_smiles,
    }

    # scs = sorted(m.sidechains, key=lambda x: len(x), reverse=True)
    scs = m.sidechains
    scds = []
    for sc in scs:
        scu = sc.umbrella

        scd = {
            'ang_vp': sc.angle_vp,
            'natoms': len(sc),
            'volume': sc.volume,
            'canonical_smiles': sc.canonical_smiles,
            'maxrank': sc.maxrank,
            "branch_rank": sc.branch_msite_rank,
            "is_hydrogen": sc.ishydrogen,
            "has_ring": sc.hasring,
            "scid": sc.scid,
            "descriptors": sc.shape_descriptors,
            "similarity_v_tips": sc.fp_similarity(tips_sc),
        }
        if scu != None:
            scud = {
                'volume': scu.volume,
                'canonical_smiles': scu.canonical_smiles,
                "descriptors": scu.shape_descriptors,
                "similarity_v_tipsu": scu.fp_similarity(tips_sc_u),
            }
            scd['scu'] = scud
        else:
            scd['scu'] = None
        scds.append(scd)

    mol_d = {
        'canonical_smiles': m.canonical_smiles,
        'volume': m.volume,
        'bone': bd,
        'sidechains': scds,
    }
    core_data = {
        'cifname': name,
        'packing_data': packingdata,
        'mol_data': mol_d
    }
    return name, core_data

tips_sc = get_tips()

for cif in sorted(glob.glob('../../cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
    if int(cif.split('/')[-1][:-4]) >= 196:
        print('working on {}'.format(cif))
        name, core_data = parser(cif, tips_sc)
        with open('{}.json'.format(name), 'w') as fp:
            json.dump(core_data, fp)

