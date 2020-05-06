import hashlib
from abc import ABCMeta
from abc import abstractmethod
from collections import OrderedDict

import anytree
from monty.json import MSONable


def hashstring(string, method=hashlib.sha256):
    return method(string.encode('utf-8')).hexdigest()


"""
Workflow
    Workflow is a python module based on fireworks, it is highly encouraged to save a snapshot of this module
    before production. (in case you changed it afterwards and there is no code for you to debug)

ParamSet
    ParamSet is a collection and 
    its document contains the arguments for the input generation (e.g. external calculators) functions of Workflow.
    e.g. ParamSet.vasp.blabla['encut']=500 eV

OCELOT DB
    this is the mongodb used to store data, hosted on lcc vm right now

Considering the outputs from OCELOT workflow can be modularized as a set of related Structure-Results pairs,
in OCELOT DB, we use the following (abstract) classes to represent the computational results:

- DbSchema
    this is an abstract class, the users should write their own subclasses to store outputs from Workflow
    this is a class representing **one** molecular/crystal structure.
        e.g. a pymatgen structure, a configuration, a rdkit molecule, or even a MolGraph. 
    this is where a set of calculations starts, that is, 
        the structural information it contains should be complete enough to be used in input generation
        as long as the ParamSet is known
    this can be used in different workflows.
    
    it must have the following attributes defined:
    self.results: OrderedDict
        this is a placeholder for results, as the workflow progress it will be filled with OcelotResult, default is an empty OrderedDict
    self.origins: [_id field of a DataScheme]
        so we can trace back to the very original input, e.g. CuratedData
        note it is a list as multiple origins are allowed, e.g. polymorphism
    self.structure:
        a structural obj from pymatgen, ocelot, rdkit, etc
    self._id:
        the hash of self.structure
        notice in MongoDB it is possible for 2 documents to have the same _id from 2 different collections
    
    we do not include self.children or self.parent here as one DbSchema can be in multiple trees, e.g. TIPGe MolGraph can
    be in TIPGe-BW tree or TIPGe-SS tree. To avoid this we have to define another hash method so 2 DbMolGraph objects can 
    have different _id but the same structural_hash.
 
- DbSchemaTree       
    this is a tree object from the anytree package, the nodes are the `_id`s of DbSchema objects
    queries can be made to this collection so a set of related DbSchema objects can be retrieved
    
    self.origin: _id field of a DataScheme
        each tree has one and only one origin
    self.workflow: str
        the name of the workflow used to 'grow' this tree
    self.paramset: str
        the name of the paramset used to 'grow' this tree
    self._id:
        hash(self.origin, self.workflow, self.paramset)
    self.root: DbSchema._id
        the root node of a tree, each node is a DbSchema._id
    self.nodes: [DbSchema._id]
        all nodes in a list, generated from one of the anytree.iterators

- OcelotResult
    designed for Web UI and ML usage, an abstract class
    it should have a relatively flat structure or at least a method to flat itself
    
    self.belongs_to: DbSchema._id
    self.workflow: str
    self.paramset: str
    self.corehours: float
    self.check(): -> bool
    self.visual: [dict]
        used in Web UI
        self.visual[0] = {'header': 'phonon bands', 'method': 'interactive_image', 'data': <data to be plot>}
    self.descriptors: dict
        a flat dict for all descriptors
    self.timestamp: str
    self._id: hash(self.timestamp)  # or just timestamp/mongo default? we want to make each OcelotResult unique


The first step of the entire Workflow is to generate the root node of DbSchema and create a DbSchemaTree, 
along with other Dbschema objects if any

    Origin + Workflow + ParamSet 
        --bootstrap--> DbSchemaTree.root + DbSchemaTree + other DbSchema (if any)
        
as the workflow progresses, DbSchemaTree may have more nodes and DbSchema.results is filled

    Workflow --run--> raw outputs on HPC + OcelotResult objects (to be filled to DbSchema)
"""


class DbSchema(MSONable, metaclass=ABCMeta):
    def __init__(self, structure, results: OrderedDict = None, origins: [str or bytes] = None, ):
        self.structure = structure
        if results is None:
            self.results = OrderedDict()
        else:
            self.results = results
        self.origins = origins

    @abstractmethod
    @property
    def _id(self) -> str or bytes:
        """hash self.structure, will this be picked up by MSONable.as_dict()?"""
        pass

    @abstractmethod
    def mongo_insert(self, *args):
        """ when inserting first check duplicates, if already existed then append origin to self.origins"""
        pass


class DbSchemaTree(MSONable):

    def __init__(self, root: anytree.node.node, origin: str or bytes, workflow: str = None, paramset: str = None):
        self.root = root
        self.origin = origin
        self.workflow = workflow
        self.paramset = paramset

    @property
    def _id(self):
        return hashstring(self.workflow + self.paramset + self.origin)

    @property
    def nodes(self):
        return anytree.iterators.preorderiter.PreOrderIter(self.root)

    @classmethod
    def from_root(cls, *args):
        # TODO create cls from a root node
        pass

    def add_child(self, parent):
        # TODO add child to a parent node
        pass


class OcelotResult(MSONable, metaclass=ABCMeta):
    def __init__(self, belongs_to, workflow, paramset, corehours, timestamp):
        self.belongs_to = belongs_to
        self.workflow = workflow
        self.paramset = paramset
        self.corehours = corehours
        self.timestamp = timestamp

    @abstractmethod
    @property
    def check(self):
        pass

    @abstractmethod
    def visual(self):
        pass

    @abstractmethod
    @property
    def descriptors(self):
        pass

    @property
    def _id(self):
        # TODO or just timestamp/mongo default? we want to make each OcelotResult unique
        return hashstring(self.timestamp)
