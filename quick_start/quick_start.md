
## Installation
0. make sure `conda` works properly
1. clone the latest stable (e.g. `v0.2`) by 
    ```bash
    git clone --single-branch --branch v0.2 git@github.com:caer200/ocelot_api.git
    ```
2. this yields a folder called `ocelot_api`, create a new venv with
    ```bash
    conda env create -f venv/env.yml
    ``` 
   or
   ```bash
   conda create --name ocelot venv/spec-file.txt
   ```
3. run `python setup.py install` there

## schema

A `MolGraph` is a graph consists of integer nodes (nodename) and edges connecting them.
For each node, there is an string attribute denotes the element.

A `MolGraph` can be partitioned into a set of `FragmentGraph`. A `FragmentGraph` contains
the information of `joints` at which fragmentation happened. The `FragmentGraph` is
used to represent different functional fragments commonly seen in functional organic
molecules.

One level above the `MolGrap` is the molecualr graph that contains
 details of bonds and basic information about the electronic system.
This is the "molecule" drawn by chemists and 
can be nicely described by the molecule class in `rdkit`. It should be noticed that
converting `MolGraph` to `rdkit.mol` is not trivial. 
We use the method from [xyz2mol](https://github.com/jensengroup/xyz2mol) by 
[Jensen Group](https://github.com/jensengroup).

With conformational information, a `MolGraph`/`FragmentGraph`
 becomes a `MolConformer`/`FragConformer` that can
be uniquely defined by the Cartesian coordinates (xyz) of its atoms. Basically, 
they are `pymatgen.Molecule` except for each `site`, there is a property `siteid`
that can be mapped to the nodes in `MoleGraph`. 

Adding periodicity to a `MolConformer` yields a `Config`, which is just a 
`pymatgen.strcuture` with no disordered sites. Disorder should be represented by
a set of weighted `Config`, as this API is used primarily for 
handling organic molecular crystals in which the number of possible configurations is limited
by molecular structure.

## Disorder



For more info see [disorder_test](../tests/disorder_test)



## Quick Start

Let's say you have a CIF file as `tipgebw.cif`. It looks like this in `Jmol`.

![tipgebw][tipgebw_jmol]

[tipgebw_jmol]: ./tipgebw.png

Such structure is not ready for any calculations as it contains disorder (at TIPS),
if you look into the cif file you will find 2 disorder groups. So the first thing we want is to clean up disorder.
```python
from ocelot.routines.pbc import CIFparser

with open('tipgebw.cif', 'r') as f:  # read cif file into string
    fs = f.read()
    
cp = CIFparser.from_cifstring(fs)  # create a parser object from cif string

clean_cif_strings = cp.get_clean_cifs_stringlist()  # a list of cif strings without disorder

with open('tipgebw_clean_0.cif', 'w') as f:  # write the first cleaned cif string into file
    f.write(clean_cif_strings[0])
```
Here we created a `CIFparser` object, then used its `get_clean_cifs_stringlist` method to 
split a cif file with 2 disorder groups into 2 cif file strings with no disorder.

Once we have a 'cleaned' cif file, we can now initiate a `Config` object to represent
the periodic structure. `Config` is a component of `ocelot` schema.
```python
from ocelot.schema.config import Config

# init Config object from cleaned cif string
tipge_config = Config.from_cifstring(clean_cif_strings[0])  

# we can remove pbc and focus on the organic molecule in this configuration
# OMol -- organic molecule -- is a list of MSites
omol_0 = tipge_config.omols[0]  
print(omol_0) 
# OMol:
# MSite: Ge (1.8797, 1.8744, 14.3509) siteid: 0
# MSite: C (2.7783, 2.5445, 12.7950) siteid: 1
# MSite: C (0.0097, 2.4446, 14.2340) siteid: 2
# ...

# most of the time we just care conjugate backbone, is also a list of MSites
bone = omol_0.backbone  
print(bone)
# Backbone:
# MSite: C (3.8686, 3.3564, 10.5508) siteid: 15
# MSite: C (4.5668, 2.4190, 9.7409) siteid: 34
# MSite: C (3.7152, 4.7050, 10.1166) siteid: 35
# ...

# we can also check sidechains, notice even a single H is considered as one sidechain
# again, sidechain is a list of MSites
sidechains = omol_0.sidechains  
print(len([sc for sc in sidechains if not sc.ishydrogen]))  # how many non-H side chains?
# 2
```
From above examples you may have noticed the structure of `ocelot` schema, 
it can be represented as something like this:
```
Config
    - OMol              |  - Ring
        - Backbone      |       - Bond
        - Sidechain     |
    |--       a list of MSites       --|
```

