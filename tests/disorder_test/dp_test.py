from ocelot.routines.disparser import DisParser
import glob
from pprint import pprint

# nodis-0 AMOZOT.cif
# nodis-1 ASEZEE.cif
# AGUHUG.cif(we cannot deal with this type)
if __name__ == '__main__':
    for i in glob.glob('x17059.cif'):  # dis-alpha
    # for i in glob.glob('ALOVOO.cif'):  # dis-gamma
    # for i in glob.glob('AGUHUG.cif'):  # dis-gamma, cannot handle
    # for i in glob.glob('AMOZOT.cif'):  # noids-0, cannot handle
    # for i in glob.glob('ASEZEE.cif'):  # dis-beta, cannot handle
        dp = DisParser.from_ciffile(i)
        r = dp.to_configs(vanilla=False, write_files=True)
        # r = dp.to_configs(vanilla=True, write_files=True)


