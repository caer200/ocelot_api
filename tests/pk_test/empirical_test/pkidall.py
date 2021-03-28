from ocelot.schema.configuration import Config
import glob
from pprint import pprint
from ocelot.task.pkidopt import pkid_ciffile, PkidError

"""
write pk params to 'packing.txt' for all cifs in ./cifs/
'./japacking.txt' shows the packing pattern identified by JA manually
"""


import pandas as pd
df_human = pd.read_csv("japacking.txt", sep="\t", header=None)
df_human.columns = ["cif", "pkid_human", "backbone", "codename"]

folder = "./cifs"

packings = []
for cif in df_human["cif"]:
    print("---")
    print("working on:", cif)
    filename = "{}/{}".format(folder, cif)
    try:
        data = pkid_ciffile(filename)
    except Exception as e:
        print("error in pkid!")
        print(str(e))
        packings.append("")
        continue
    packing = []
    for k in data.keys():
        packing.append(data[k]["packing"])
    packings.append("__".join(sorted(set(packing))))
df_human["pkid_ocelot"] = packings
df_human.to_csv("pkid.csv")



# f = open('packing.txt', 'w')
# f.write(
#     '# cifname packing  n_close_azm_and_parallel  n_close_azm_and_notparallel  n_close_vertical_and_parallel  n_close_vertical_and_notparallel  n_parallel_and_overlap  n_notparallel_and_overlap  \n')
# for cif in sorted(glob.glob('./cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
#     print('working on: ' + cif)
#     n, packingdata, lbone = inspect_cif(cif)
#     packingdata_list = [p for p in packingdata.values()]
#     outstring = ''
#     for refp in packingdata_list:
#         packingstring = '-'.join([str(va) for va in refp.values()])
#         outstring += packingstring + '@'
#     f.write("{}\t\t{}\t\t{:.3f} \n".format(n, outstring, lbone))
# f.close()
