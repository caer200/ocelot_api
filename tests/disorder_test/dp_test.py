from ocelot.routines.disparser import DisParser
import glob

# nodis-0 AMOZOT.cif
# nodis-1 ASEZEE.cif
# dis-0 x17059.cif
# dis-1 tipge.cif
# dis-2 AGUHUG.cif(we cannot deal with this type) ALOVOO.cif
if __name__ == '__main__':
    # dp = DisParser.from_ciffile('gebw.cif')
    # for i in glob.glob('x17059.cif'):
    for i in glob.glob('ALOVOO.cif'):
        dp = DisParser.from_ciffile(i)
        results = dp.to_configs(write_files=True)
    # pprint(results)
    #     a = dp.classify()
    #     print(a, i)


