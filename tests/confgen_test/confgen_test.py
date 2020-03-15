from ocelot.task.confgen import *

smiles = {
    'biphenyl': 'c1ccccc1-c2ccccc2',
    'dodecane': 'CCCCCCCCCCCC',
    'harmine': 'COc1ccc2c(c1)[nH]c3c(C)nccc23',
    'coronene': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67',
}

for method in ['cluster', 'prune']:
    for n, smi in smiles.items():
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        cg = ConfGen(m, n, prune_in_gen=0.5, rngseed=42, nconf=20)
        m_after, e_after = cg.genconfs(flavor=method, write_confs=True, prune_max_conformer=20, rmsd_threshold=0.5)
        print('working on {}: {}'.format(method, n), len(e_after))

# TODO there seems to be a bug with biphenyl, it gives 2 confs with rmsd ~0.9, this seems to be related
#  to the trapped minimization of rmsd which was induced by atom index mapping.
#  looping over all mappings seems unreasonable as mol symm gets higher
#  we can try e.g. similarity measure based on atomic environment to get rid of these cases
