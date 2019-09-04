"""
this contains the schema for dealing with conjugate organics, both molecules and beyond

hierarchy as:
element -- > msite --> msitelist --> bond --> ring/sidechain/backbone --> omol --> dimer/config

one should try not to init anything lower than omol individually as they rely heavily on the omol (site_id)
"""
