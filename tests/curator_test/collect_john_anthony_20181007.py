import os
import glob

from ocelot.curator.Contribution import Contribution
from ocelot.curator.DataSchema import *

ACCESS_PROVIDER = JohnAnthony_UKY
ACCESS_DATE = datetime.date(2018, 10, 7)
ACCESS_LIC = None
data_access = DataAccess(ACCESS_PROVIDER, ACCESS_DATE, ACCESS_LIC)

dir_path = os.path.dirname(os.path.realpath(__file__))

CONTRIBUTION_FOLDER_BASENAME = 'John_Anthony_20181007'

CONTRIBUTION_FOLDER = '/home/ai/tmp/community_json/community_cifs/contributions/John_Anthony_20181007'
COLLECT_FOLDER = '{}/collect_{}'.format(dir_path, CONTRIBUTION_FOLDER_BASENAME)
CURATE_FOLDER = '{}/curate_{}'.format(dir_path, CONTRIBUTION_FOLDER_BASENAME)

john_anthony_contribution_20181007 = Contribution(
    data_access, CONTRIBUTION_FOLDER
)

def john_anthony_vis_pkid(abspath):
    common_path = os.path.commonpath([CONTRIBUTION_FOLDER, abspath])
    relative_path = os.path.relpath(abspath, common_path)
    d = OrderedDict()
    items = relative_path.split('/')
    if len(items) == 3:
        pk, chrom, basename = items
    else:
        pk, basename = items
        chrom = 'unknown'
    d['visually_inspected_packing'] = pk
    d['visually_inspected_chromophore'] = chrom
    return d

if __name__ == '__main__':
    """
    ./contributions/John_Anthony_20181007/2D_brickwork/Anthradithiophene/x14433-not-good_2ch2oh_tips_adt.cif
        a 4x1x1 supercell was saved as a unitcell, there is only one oxygen atom per molecule, no hydrogen was fitted
        this cif file will not be collected
    """
    MLAB_URI = "mongodb://qai222:caer200@ds047782.mlab.com:47782/ocelot_qai?retryWrites=false"
    LCC_VM_URI = "mongodb://ocelot:caer200@10.33.28.79:27017/"
    # import time
    # ts1 = time.time()
    # collected_filepaths, raw_data_list, jsonfilepaths = john_anthony_contribution_20181007.collect_with_log(dir_path, COLLECT_FOLDER, john_anthony_vis_pkid)
    # ts2 = time.time()
    # print('$$$ collection step took {}'.format(ts2 - ts1))  # 19

    # jsonfilepaths = glob.glob('{}/*.json'.format(COLLECT_FOLDER))
    # raw_data_list = [DataEntry.from_jsonfile(jf) for jf in jsonfilepaths]
    # for rd in raw_data_list:
    #     rd.mongo_insert(MLAB_URI, 'ocelot_qai')

    # import time
    # ts1 = time.time()
    # curated_data_list, classes = john_anthony_contribution_20181007.curate_all(dir_path, raw_data_list, CURATE_FOLDER)
    # ts2 = time.time()
    # print('$$$ curate step took {}'.format(ts2 - ts1))  # 1331

    # push to mongo
    c_jsons = sorted(list(glob.glob('{}/*/*.json'.format(CURATE_FOLDER))))
    c_data_list = [CuratedData.from_jsonfile(jf) for jf in c_jsons]

    for cd in c_data_list:
        cd.mongo_insert(MLAB_URI, 'ocelot_qai')
