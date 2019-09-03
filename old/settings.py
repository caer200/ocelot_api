# global variables
SYSPATH = dict(
    # a folder contains all cif files need to be calculated
    CIF_REPO='/scratch/qai222/proj/m18209/raw',
    # for vasp PAW
    PAW_PATH='/scratch/qai222/database/PAW_5.4.1/',
    # a folder for compiled binaries
    BIN_PATH='/scratch/qai222/database/aflow/',
    # a folder for all runs
    RUN_PATH='/scratch/qai222/proj/m18209/run',
    ANA_PATH='/scratch/qai222/proj/m18209/analysis',
    ANA_TEMP='/scratch/qai222/proj/m18209/',
)



# TD0 = {'what':'analyze raw cif file', 'dependents':[], 'bin':'cifparser', 'number':0}
# TD1 = {'what':'single mol relaxation', 'dependents':[0], 'bin':'vasp', 'number':1}
# TD2 = {'what':'fixbox relaxation', 'dependents':[0], 'bin':'vasp', 'number':2}
# TD3 = {'what':'fixbox scf', 'dependents':[2], 'bin':'vasp', 'number':3}
# TD4 = {'what':'fixbox bands', 'dependents':[0, 3], 'bin':'vasp', 'number':4}
# TD5 = {'what':'no constraint relax', 'dependents':[0], 'bin':'vasp', 'number':5}
# TD6 = {'what':'no constraint scf', 'dependents':[5], 'bin':'vasp', 'number':6}
# TD7 = {'what':'no constraint bands', 'dependents':[0, 6], 'bin':'vasp', 'number':7}
# TD8 = {'what':'dimer coupling', 'dependents':[0], 'bin':'gauss', 'number':8}
# TD9 = {'what':'nca opt', 'dependents':[0], 'bin':'gauss', 'number':9}
# TD10 = {'what':'reorganization energy', 'dependents':[9], 'bin':'gauss', 'number':10}
# TD11 = {'what':'uv vis singlet', 'dependents':[10], 'bin':'gauss', 'number':11}
#
#
#

