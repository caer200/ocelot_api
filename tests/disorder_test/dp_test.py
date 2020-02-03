from ocelot.routines.disparser import DisParser

# TODO deal with cif file from jmol, that is, no _atomic_site_label field
if __name__ == '__main__':
    from pprint import pprint
    dp = DisParser.from_ciffile('gebw.cif')
    results = dp.to_configs(write_files=True)
    pprint(results)
