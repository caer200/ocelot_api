import json
import re
import warnings
from pprint import pprint
from pymatgen.io.cif import CifFile

with open('disorder.json', 'r') as f:
    cifs = json.load(f)
keys = list(cifs.keys())

def get_pmg_dict(cifstring):
    cifdata = CifFile.from_string(cifstring).data
    idnetifiers = list(cifdata.keys())
    if len(idnetifiers) > 1:
        warnings.warn('W: find more than 1 structures in this cif file!')
    elif len(idnetifiers) == 0:
        warnings.warn('W: no structure found by pymatgen parser!')
    identifier = idnetifiers[0]
    pymatgen_dict = list(cifdata.items())[0][1].data
    return identifier, pymatgen_dict

def check_label_disorder(pmg_dict):
    """
    check if any lable is anything other than ^\w+$

    :param pmg_dict:
    :return:
    """
    is_label_disorder = False
    for label in pmg_dict['_atom_site_label']:
        if not re.search(r"^\w+$", label):
            is_label_disorder = True
            break
    return is_label_disorder


def check_site_consistency(pmg_dict):
    labels = pmg_dict['_atom_site_label']
    elements = [re.match(r"^[a-zA-Z]+", label).group(0) for label in labels]






# pprint(pymatgen_dict['_atom_site_label'])
# for k in keys:
#     cifstring = cifs[k]
#     csdid, pdict = get_pmg_dict(cifstring)
#     isdis = check_label_disorder(pdict)
#     if not isdis:
#         print(csdid)
#         with open('{}.cif'.format(csdid), 'w') as f:
#             f.write(cifstring)
