from pprint import pprint

from ocelot.curator.Contribution import *

MLAB_URI = "mongodb://qai222:caer200@ds047782.mlab.com:47782/ocelot_qai?retryWrites=false"

rawdata = RawData.from_mongo(
    mongo_query='John_Anthony-UoK-2018-10-07--k01074',  # if str is provided it searches _id
    mongo_uri=MLAB_URI,
    database_name='ocelot_qai'
)
pprint(rawdata.data_properties)

curateddata = CuratedData.from_mongo(
    mongo_query={"data_properties.rawdataid": "John_Anthony-UoK-2018-10-07--k01074"},
    mongo_uri=MLAB_URI,
    database_name='ocelot_qai',
)

c: Config = curateddata.data_content['configuration']
mc: MolConformer = c.molconformers[0]
print(type(mc.backbone.pfit_vp))  # <class 'numpy.ndarray'>
