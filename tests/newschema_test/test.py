from ocelot.schema.conformer.ocelotsite import OcelotSite
from ocelot.schema.conformer.ocelotsitelist import OcelotSiteList


osl = OcelotSiteList.from_file('tipgpn.xyz')
rd, smiles = osl.to_rdmol(0, expliciths=False)
print(smiles)


# s1 = OcelotSite('Se')
# ps = s1.to_pymatgen_site()
# s2 = OcelotSite.from_pymatgen_site(ps)

# hash(s1)
# print(s1)
# print(ps)
# print(s2)
#
#
# osl = OcelotSiteList([s1, s2])
# print(set([s1, s2]))
# print(len(osl))
# print(s1 == s2)
# print(s2 in osl)
# print(osl.idxlist)

