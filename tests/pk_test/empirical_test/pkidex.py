from ocelot.task.pkidopt import pkid_ciffile
from pprint import pprint
"""
write pk params to 'packing.txt' for all cifs in ./cifs/
'./japacking.txt' shows the packing pattern identified by JA manually
"""

import pandas as pd

df_human = pd.read_csv("japacking.txt", sep="\t", header=None)
df_human.columns = ["cif", "pkid_human", "backbone", "codename"]

folder = "./cifs"

packings = []
exceptions = [
    # 31,  # algo right, human right
    # 51,  # human wrong
    # 52,  # human wrong
    # 68,  # human wrong
    # 113,  # human wrong
    # 122,  # human wrong
    # 126,  # human wrong
    # 137,  # algo wrong due to severe disorder leading to wrong bone.cif
    # 142,  # human wrong
    # 164,  # human wrong
    # 167,  # algo right, human right
    # 170,  # algo right, human right
    # 174,  # algo right, human right
    # 175,  # algo right, human right
    # 176,  # algo right, human right  #TODO the "ribbon" label from J.A. seems to be slipped-stack + 1D structure
    # 180,  # human wrong
    # 181,  # algo right, human right
    # 182,  # algo right, human right
    # 184,  # human wrong
    # 185,  # human wrong
    # 187,  # algo right, human right
    # 189,  # algo right, human right
    # 191,  # algo right, human right
    # 193,  # algo right, human right
    # 195,  # algo right, human right
]
for r in df_human.to_dict("records"):
    cif = r["cif"]
    cifid = int(cif[:-4])
    if cifid not in exceptions:
        continue
    print("---")
    print("J.A. believes this is: {}".format(r["pkid_human"]))
    print("working on:", cif)
    filename = "{}/{}".format(folder, cif)
    data = pkid_ciffile(filename)
    pprint(data)
    break
