import calendar
import glob
import os
import re
import zipfile

from docx2csv import extract_tables

from ocelot.curator.Contribution import Contribution
from ocelot.curator.Contribution import createdir
from ocelot.curator.DataSchema import *

ACCESS_PROVIDER = MikeHaley_UO
ACCESS_DATE = datetime.date(2018, 4, 23)
ACCESS_LIC = None
data_access = DataAccess(ACCESS_PROVIDER, ACCESS_DATE, ACCESS_LIC)

dir_path = os.path.dirname(os.path.realpath(__file__))

CONTRIBUTION_FOLDER_BASENAME = 'Mike_Haley_20180423'

CONTRIBUTION_FOLDER = '/home/ai/tmp/community_json/community_cifs/contributions/Mike_Haley_20180423'
COLLECT_FOLDER = '{}/collect_{}'.format(dir_path, CONTRIBUTION_FOLDER_BASENAME)
CURATE_FOLDER = '{}/curate_{}'.format(dir_path, CONTRIBUTION_FOLDER_BASENAME)

mike_haley_contribution_20180423 = Contribution(data_access, CONTRIBUTION_FOLDER)

DOCXFILE = CONTRIBUTION_FOLDER + '/Haley_IF_Structures.docx'


def getts():
    return calendar.timegm(time.gmtime())


def mike_haley_report(docxfile=DOCXFILE, wdir=dir_path + '/scratch-{}/'.format(getts())):
    """
    mike provided a docxfile which is just a table of 3 columns:
    emf image of molecule, crystal structure codename, published paper (prior to 20180423)
    """
    createdir(wdir)
    whereami = os.getcwd()
    os.chdir(wdir)
    z = zipfile.ZipFile(docxfile)
    all_files = z.namelist()
    images = filter(lambda x: x.startswith('word/media/'), all_files)
    for imgname in images:
        img = z.open(imgname).read()
        f = open(os.path.basename(imgname), 'wb')
        f.write(img)
        f.close()
    z.close()
    emffiles = sorted(list(glob.glob('{}/*.emf'.format(wdir))),
                      key=lambda x: int(os.path.basename(x).strip('.emf').strip('image')))
    img_column = []
    i = 1
    for emffile in emffiles:
        print('converting {}/{}'.format(i, len(emffiles)))
        i += 1
        cmd = 'libreoffice --headless --convert-to eps {}; sleep 2'.format(emffile)
        os.system(cmd)
        # the sleep has to be here... see
        # https://bugs.launchpad.net/ubuntu/+source/libreoffice/+bug/1777285
        epsfile_base = os.path.basename(emffile).strip('.emf') + '.eps'
        epsfile = '{}/{}'.format(wdir, epsfile_base)
        with open(epsfile, 'r') as f:
            b = f.read()
            img_column.append(b)

    print(len(img_column))
    table = extract_tables(docxfile)[0][2:]
    summary = OrderedDict()
    print(len(table))
    for i in range(len(table)):
        _, mh_codename, publication = [item.decode() for item in table[i]]
        mh_codename = mh_codename.lower()
        if 'unpublish' in publication:
            publication = None
        if len(mh_codename.split()) > 1:
            for codename in mh_codename.split():
                summary[codename] = [img_column[i], publication]
        else:
            summary[mh_codename] = [img_column[i], publication]
    summary_jsonpath = '{}/{}-summary.json'.format(wdir, data_access)
    with open(summary_jsonpath, 'w') as f:
        json.dump(summary, f)
    os.chdir(whereami)
    return summary, summary_jsonpath


if __name__ == '__main__':
    # summary, summary_jsonpath = mike_haley_report()
    summary_jsonpath = '/home/ai/ocelot_api/tests/curator_test/scratch-1588142735/Mike_Haley-UoO-2018-04-23-summary.json'
    with open(summary_jsonpath, 'r') as f:
        summary = json.load(f)


    def assign_properies(abspath):
        codename = os.path.basename(abspath).strip('.cif').lower()
        codename = re.match(r"[a-z]+\d+", codename).group(0)
        d = OrderedDict()
        try:
            props = summary[codename]
        except KeyError:
            return d
        d['mike_haley_molecule_structure_plot_eps'] = props[0]
        d['publication_information_upto_20180423'] = props[1]
        return d


    MLAB_URI = "mongodb://qai222:caer200@ds047782.mlab.com:47782/ocelot_qai?retryWrites=false"

    # LCC_VM_URI = "mongodb://ocelot:caer200@10.33.28.79:27017/"
    # ts1 = time.time()
    # collected_filepaths, raw_data_list, jsonfilepaths = mike_haley_contribution_20180423.collect_with_log(dir_path, COLLECT_FOLDER, assign_properies)
    # ts2 = time.time()
    # print('$$$ collection step took {}'.format(ts2 - ts1))  # 19

    jsonfilepaths = glob.glob('{}/*.json'.format(COLLECT_FOLDER))
    raw_data_list = [DataEntry.from_jsonfile(jf) for jf in jsonfilepaths]
    for rd in raw_data_list:
        rd.mongo_insert(MLAB_URI, 'ocelot_qai')

    import time

    ts1 = time.time()
    c_data_list, classes = mike_haley_contribution_20180423.curate_all(dir_path, raw_data_list, CURATE_FOLDER)
    ts2 = time.time()
    print('$$$ curate step took {}'.format(ts2 - ts1))  # 1331
    #
    #     # push to mongo
    #     c_jsons = sorted(list(glob.glob('{}/*/*.json'.format(CURATE_FOLDER))))
    #     c_data_list = [CuratedData.from_jsonfile(jf) for jf in c_jsons]
    #
    for cd in c_data_list:
        cd.mongo_insert(MLAB_URI, 'ocelot_qai')
