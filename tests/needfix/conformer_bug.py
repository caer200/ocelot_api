"""
mol conformer to_rdmol failed for csd_TIQBUY
it seems there are more than one bugs here:
    - the xyz from molconformer has two molecules
    - the mc.to_rdmol failed (may) due to mem explosion in
        `valences_list = list(itertools.product(*valences_list_of_lists))`
"""