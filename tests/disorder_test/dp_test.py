from ocelot.routines.disparser import DisParser
import glob
from pprint import pprint

# nodis-0 AMOZOT.cif
# nodis-1 ASEZEE.cif
# dis-2 AGUHUG.cif(we cannot deal with this type)
if __name__ == '__main__':
    # dp = DisParser.from_ciffile('gebw.cif')
    for i in glob.glob('x17059.cif'):  # dis-0
    # for i in glob.glob('tipge.cif'):    # mis-classified to dis-1, should be dis-0
    # for i in glob.glob('ALOVOO.cif'):  # dis-2
    # for i in glob.glob('tipge_conf0.cif'):    # nodis-0
        dp = DisParser.from_ciffile(i)
        dp.to_configs(write_files=True)


