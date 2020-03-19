"""
Crystal:
    '_id':
        source + identifier  # e.g. 'csd_ASFNC'
    'hash':
        elements + lattice_matrix  # e.g. hash('c', 'h', 9.2, 7.2, 17.5, ...)
    'identifier':
        string 
    'source':
        string  
    'cif_string': 
        string  
    'elements': 
        (str)  # e.g. ('c', 'h', ...), must be sorted
    'disorder_location':
        'bone' or 'sidechain' or 'not sure' or 'no disorder
    'disorder_class':
        string  # from DisParser
    'disordered_structure_wrapped': 
        pmg.structure.as_dict()
    'configurations':
        [Configuration._id]
    'major_configuration':
        Configuration._id

    --- universal start---
    'lattice': 
        pmg.lattice.as_dict()
    'smiles':
        string
    'hsmiles':
        string
    'chromophore':
        Chromophore._id
    --- universal end  ---

    'fin_data':
        a dictionary for storing calculated data, keys can be 
            "packing_data", "hop_data", "ebs", "lem", "tem", "cohesive", ...

Configuration:
    '_id':
        id(CrystalEntry) + config_id  # e.g. 'csd_ASFNC_conf0'
    'hash':
        elements + lattice_matrix + occupancy + config_id
    'calculation_stage':
        string  # denotes how many calculations have been done for this config
    'unwrap_structure':
        pmg.structure.as_dict()  # this is the unwrap_clean_pstructure
    'molconformers':
        [MolConformer._id]  # use rmsd to get rid of identical conformers, see BasicConformer.compare()
    'z':
        int
    'occupancy':
        float
    'backbone_structure':
        pmg.structure.as_dict()  # this is the clean bone structure with hydrogen terminations
    'Crystal':
        Crystal._id
    'config_id':
        int
    'packing_data':
        dict  # returned by packing identification
    'hop_data':  # can be geodist, zindo, dft
        json returned by hop.hopping
    'fb_opt_structure':
        pmg.structure.as_dict()
    'ar_opt_structure':
        pmg.structure.as_dict()
    'fb_opt_ebs':
        dict  # keys are eigens, line em, tem and something to store the graph
    'ar_opt_ebs':
        dict  # keys are eigens, line em, tem and something to store the graph
    'energetics':
        dict  # keys are fb_opt, ar_opt, mol_in_box, cohesive

MolConformer:
    '_id':
        id(Configuration) + mc_index
    'mc_index':
        int
    'hash':
        id(MolConformer)
    'Configuration':
        Configuration._id
    'GeoBoneConf':
        BackboneConformer.as_dict()
    'GeoBoneConf_descriptors':
        dict
    'GeoScsConf':
        [SidechainConformer.as_dict()]
    'GeoScsConf_descriptors':
        [dict]
    'ChromBoneConf':
        BackboneConformer.as_dict()
    'ChromBoneConf_descriptors':
        dict
    'ChromScsConf':
        [SidechainConformer.as_dict()]
    'ChromScsConf_descriptors':
        [dict]

ChromophoreConf:
    '_id':
        hsmiles + index
    'index':
        int  # give by confgen
    'geo':  # original generated geometry
        pmg.Molecule.as_dict()
    'anion_geo':
        pmg.Molecule.as_dict()
    'cation_geo':
        pmg.Molecule.as_dict()
    'AIP':
        float
    'AEA':
        float
    'reorg':
        float
    'tddft':
        dict  # transitions
    's0s1':
        float
    's0t1':
        float
    's0t1deltascf':
        float
    'Chromophore':
        Chromophore._id

Chromophore:
    '_id':
        hsmiles  # or other method hashing a molecule
    'hash':
        id(Chromophore)
    'ChromophoreConfs':
        [ChromophoreConf._id]
"""
