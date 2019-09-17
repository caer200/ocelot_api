
## Installation
1. clone the release branch (e.g. `v0.01`) by 
    ```bash
    git clone --single-branch --branch v0.01 git@github.com:caer200/ocelot_api.git
    ```
2. create a virtual env and install
   ```bash
   conda create --name ocelot
   conda activate ocelot
   conda install python
   conda install --file spec-file.txt
   ```
   in case you wonder, `conda env create -f env.yml` may not work as part of `pymatgen` 
   uses pre-link scripts, see [this issue](https://github.com/materialsproject/pymatgen/issues/555).

3. run `python setup.py install`


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

