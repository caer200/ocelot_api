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

# def check_label_disorder(pmg_dict):
#     """
#     check if any lable is anything other than ^\w+$
#
#     :param pmg_dict:
#     :return:
#     """
#     is_label_disorder = False
#     for label in pmg_dict['_atom_site_label']:
#         if not re.search(r"^\w+$", label):
#             is_label_disorder = True
#             break
#     return is_label_disorder

class AtomLabel:
    def __init__(self, label:str):
        self.label = label

        tmplabel = self.label
        self.element = re.findall(r"^[a-zA-Z]+", tmplabel)[0]

        tmplabel = tmplabel.lstrip(self.element)
        self.index = re.findall(r"^\d+", tmplabel)[0]

        tmplabel = tmplabel.lstrip(self.index)
        self.third_label = tmplabel
        # self.third_label = re.findall(r"^\w+", tmplabel)
        # if len(self.third_label) == 0:
        #     self.third_label = None
        # else:
        #     self.third_label = self.third_label[0]

        # disorderlabel = re.findall(r"[a-zA-Z]$", self.label)
        # if disorderlabel:
        #     self.disorderlabel = disorderlabel[0]
        # else:
        #     self.disorderlabel = None

    def is_similar(self, other):
        return self.element == other.element and self.index == other.index

def get_thirdlabels(pmg_dict):
    labels = pmg_dict['_atom_site_label']
    als = [AtomLabel(lab) for lab in labels]
    return [al.third_label for al in als]

# def check_site_disorderlabel(pmg_dict):
#     labels = pmg_dict['_atom_site_label']
#     als = [AtomLabel(lab) for lab in labels]
#     if any(al.disorderlabel for al in als):
#         return False
#     return True
    # for i in range(len(als)):
    #     al_i = als[i]
    #     print(al_i.element, al_i.index, al_i.disorderlabel)
    #     for j in range(i+1, len(als)):
    #         al_j = als[j]
    #         if al_i.is_similar(al_j) and
    #
    #
    # for al1 in als:
    #     for al2 in als:

    # elements = [re.match(r"^[a-zA-Z]+", label).group(0) for label in labels]

for k in keys:
    cifstring = cifs[k]
    csdid, pdict = get_pmg_dict(cifstring)
    # if csdid == 'ASIXEH':
    thirdlabels = get_thirdlabels(pdict)
    if any(tl is not None for tl in thirdlabels):
        print(csdid)
        print(set([tl for tl in thirdlabels if tl is not None]))
    # isdis = check_label_disorder(pdict)
    # if not isdis:
    #     if check_label_disorder(pdict):
    #         print(csdid)
        with open('{}.cif'.format(csdid), 'w') as f:
            f.write(cifstring)
