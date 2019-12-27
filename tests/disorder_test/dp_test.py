from ocelot.routines.disparser import DisParser
if __name__ == '__main__':
    from pprint import pprint
    dp = DisParser.from_ciffile('gebw.cif')
    results = dp.to_configs(write_files=True)
    pprint(results)
