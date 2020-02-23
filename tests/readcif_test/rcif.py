from ocelot.task.readcif import *
import glob
for cif in glob.glob('tipge-*.cif'):
    print('working on {}'.format(cif))
    rc = ReadCif.from_ciffile(cif)
    print(rc.properties)
    # print(rc.configs[0])
    # data = rc.where_is_disorder()
    # print(cif)
    # for i in data.keys():
    #     print('imol = {}'.format(i), data[i].conformer_properties)
    print('done!')
