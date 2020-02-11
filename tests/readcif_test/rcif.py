from ocelot.task.readcif import *
import glob
for cif in glob.glob('tip*.cif'):
    rc = ReadCif.from_ciffile(cif)
    data = rc.where_is_disorder()
    print(cif)
    for i in data.keys():
        print('imol = {}'.format(i), data[i].conformer_properties)
