import datetime
import json
from collections import OrderedDict

import license
import pymongo
from monty.json import MSONable
from monty.json import MontyDecoder
from monty.json import MontyEncoder


class DataSchemaError(Exception): pass


def string_acronym(s):
    return "".join(e[0] for e in s.split())


class DataProvider(MSONable):
    def __init__(self, name: str, institution: str, url: str):
        """
        this can represent a research group, a journal, or a database
        """
        self.name = name
        self.institution = institution
        self.url = url

    def __repr__(self):
        name = '_'.join(self.name.split())
        institution = string_acronym(self.institution)
        return '{}-{}'.format(name, institution)

    @classmethod
    def from_publication_doi(cls, doi):
        pass


class DataAccess(MSONable):  # provider + date + lic
    def __init__(self, sharedby: DataProvider, accessdate: datetime.date, license: license.base.License or None):
        self.sharedby = sharedby
        self.accessdate = accessdate
        self.license = license

    def __repr__(self):
        dp = self.sharedby.__repr__()
        date = self.accessdate.isoformat()
        return '-'.join([dp, date])

    def as_dict(self) -> dict:
        d = OrderedDict()
        d['sharedby'] = self.sharedby.as_dict()
        d['accessdate'] = self.accessdate.isoformat()
        try:
            d['license'] = self.license.name
        except AttributeError:
            d['license'] = None
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(
            DataProvider.from_dict(d['sharedby']),
            datetime.date.fromisoformat(d['accessdate']),
            license.find(d['license'])
        )


class DataEntry(MSONable):
    def __init__(
            self,
            data_content,
            data_access: DataAccess,
            _id: str or bytes,
            data_properties=None
    ):
        self.data_content = data_content
        self.data_access = data_access
        self._id = _id
        if data_properties is None:
            data_properties = dict()
        self.data_properties = data_properties

    def mongo_insert(self, mongo_uri, database_name, collection_name=None, overwrite=False):
        """
        insert RawData as an entry into a mongodb using pymongo

        :param mongo_uri: port level is enough
        :return:
        """
        if collection_name is None:
            collection_name = self.__class__.__name__
        database = pymongo.mongo_client.MongoClient(host=mongo_uri)[database_name]

        collection = database[collection_name]
        existed = collection.find_one({'_id': self._id})
        inserted = False
        if existed is None:
            result = collection.insert_one(json.loads(self.to_json()))
            if result.acknowledged:
                inserted = True
        else:
            if overwrite:
                result = collection.replace_one({'_id': self._id}, self.as_dict())
                if result.acknowledged:
                    print('{}: previously existed {} was replaced'.format(self._id, result.matched_count))
                    inserted = True
            else:
                print('{}: previously existed and was NOT updated'.format(self._id))
        return inserted

    @classmethod
    def from_jsonfile(cls, jsonfile):
        with open(jsonfile, 'r') as prev_json:
            prev_jsons = prev_json.read()
        return json.loads(prev_jsons, cls=MontyDecoder)

    def to_jsonfile(self, jsonfile):
        with open(jsonfile, 'w') as cdata_json:
            json.dump(self.as_dict(), cdata_json, cls=MontyEncoder)
        with open(jsonfile, 'r') as cdata_json:
            s = cdata_json.read()
        if not s == self.to_json():
            raise DataSchemaError('rawjsonfile is different from to_json')
        return jsonfile

    @classmethod
    def from_json(cls, s):
        return cls.from_dict(json.loads(s))

    @classmethod
    def from_mongo(cls, mongo_query, mongo_uri, database_name, collection_name=None):
        if collection_name is None:
            collection_name = cls.__name__
        if isinstance(mongo_query, str):
            mongo_query = {'_id': mongo_query}
        database = pymongo.mongo_client.MongoClient(host=mongo_uri)[database_name]

        collection = database[collection_name]
        cursor = collection.find(mongo_query)
        if cursor.count() > 1:
            raise DataSchemaError('return multiple results w. query: {}'.format(mongo_query))
        elif cursor.count() == 0:
            raise DataSchemaError('return no results w. query: {}'.format(mongo_query))
        else:
            existed = cursor[0]
        return cls.from_dict(existed)


class CuratedData(DataEntry):

    def mongo_insert(self, mongo_uri, database_name, collection_name=None, overwrite=False,
                     rawdatacheck=True, rawdatacollection='RawData'):
        """
        check if raw data is there
        """
        if collection_name is None:
            collection_name = self.__class__.__name__
        database = pymongo.mongo_client.MongoClient(host=mongo_uri)[database_name]

        if rawdatacheck:
            rawcoll = database[rawdatacollection]
            rawid = self.data_properties['rawdataid']
            if rawcoll.find_one({'_id': rawid}) is None:
                print('this curated data relies on non-existing raw data: {}'.format(rawid))
                return False

        collection = database[collection_name]
        existed = collection.find_one({'_id': self._id})
        inserted = False
        if existed is None:
            result = collection.insert_one(json.loads(self.to_json()))
            if result.acknowledged:
                inserted = True
        else:
            if overwrite:
                result = collection.replace_one({'_id': self._id}, self.as_dict())
                if result.acknowledged:
                    print('{}: previously existed {} was replaced'.format(self._id, result.matched_count))
                    inserted = True
            else:
                print('{}: previously existed and was NOT updated'.format(self._id))
        return inserted


class RawData(DataEntry): pass


JohnAnthony_UKY = DataProvider('John Anthony', 'University of Kentucky', 'https://chem.as.uky.edu/users/anthony')
MikeHaley_UO = DataProvider('Mike Haley', 'University of Oregon', 'https://haleylab.uoregon.edu/')
CSD_CCDC = DataProvider('Cambridge Structural Database', 'Cambridge Crystallographic Data Centre',
                        'https://www.ccdc.cam.ac.uk/')
SeanParkin_UKY = DataProvider('Sean Parkin', 'University of Kentucky', 'https://xray.uky.edu/')
#
# if __name__ == '__main__':
#     ACCESS_PROVIDER = JohnAnthony_UKY
#     ACCESS_DATE = datetime.date(2018, 10, 7)
#     ACCESS_LIC = None
#     data_access = DataAccess(ACCESS_PROVIDER, ACCESS_DATE, ACCESS_LIC)
#     rd = DataEntry('lalala', data_access, 'test')
#     MLAB_URI = "mongodb://qai222:caer200@ds047782.mlab.com:47782/ocelot_qai?retryWrites=false"
#     LCC_VM_URI = "mongodb://ocelot:caer200@10.33.28.79:27017/"
#     inserted = rd.mongo_insert(MLAB_URI, 'ocelot_qai')
#     print(inserted)
#
