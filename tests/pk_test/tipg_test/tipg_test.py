from ocelot.task.pkid import pkid_ciffile as pkid
from ocelot.task.pkidopt import pkid_ciffile as pkidopt
import glob
"""
this reads cif files from tipge-*.cif, identify packing patterns and write backbone-only configs
"""


if __name__ == '__main__':

    for cif in sorted(glob.glob('./tipge-*.cif')):
        fname = cif.split('/')[-1][:-4]
        packingdata = pkid(cif)
        print("working on", cif)
        print("-> old method:")
        print("-> {} {}\n".format(fname, packingdata[0]['packing']))
        packingdata_opt = pkidopt(cif)
        print("-> opt method:")
        print("-> {} {}\n".format(fname, packingdata_opt[0]['packing']))
