# OCELOT API

API for dealing conjugate organic molecules.

Current doc is available [here](https://caer200.github.io/ocelot_api/).

### install

within a new `conda env`, do either

```
conda install --name ocelot_api --file spec-file.txt
```

or 

```
conda env create -f environment.yml
```

### doc built

use `sphinx` with
```
sphinx-apidoc . --full -o doc -H 'ocelot' -A 'Ai, Qianxiang' -V '0.01'
```

copy html to /docs if using git pages

better method: https://daler.github.io/sphinxdoc-test/includeme.html
