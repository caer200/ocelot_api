from pymongo import MongoClient
from dimers_pair import *
import json


ZINDOLIB = '/home/vinayak/ZINDO'
ZINDOBIN = '/home/vinayak/ZINDO/zindo'
ZINDOCTBIN = '/home/vinayak/ZINDO/zindo-ct'


client = MongoClient()
db = client['ocelot_ui']
coll = db['repo_entry_withowned']
cursor = coll.find()


count = 0

while count < cursor.collection.count_documents({}):
    entry = cursor[count]
    coupling_dft = entry['dimers']['hole_transport']['electronic_coupling']
    config = Config.from_cifstring(entry['crystal']['ionic_structure']['cif'])

    mols = len(config.omols)
    if not mols > 1:
        dimer_array = config.get_dimers_array(fast=True)
        dimer_unique = get_unique_dimers(dimer_array, mols)
        os.chdir('/home/vinayak/OCELOT/ocelot_api/compare2/')
        os.system('mkdir '+str(count))
        couplings = run_zindo(dimer_unique, '/home/vinayak/OCELOT/ocelot_api/compare2/'+str(count)+'/')
        homo_couplings = [count, coupling_dft, couplings]
        os.chdir('/home/vinayak/OCELOT/ocelot_api/compare2/')
        with open('compare.dat','a') as F:
            json.dump(homo_couplings,F)
            F.write('\n')

    if not cursor.alive:
        cursor = coll.find()

    count += 1
