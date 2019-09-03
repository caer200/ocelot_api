import json
import glob
import pandas as pd
from pprint import pprint

def parse_json(fn):
    with open(fn, 'r') as f:
        d = json.load(f)
    return d

def emprical_packing(packingdata):
    packing = [packingdata[k]['packing'] for k in packingdata.keys()]
    if 'brickwork' in packing:
        p = 'brickwork'
        c = 'red'
    elif 'slipped_stack' in packing:
        p = 'slipped_stack'
        c = 'blue'
    elif 'herringbone' in packing:
        p = 'herringbone'
        c = 'gold'
    elif 'dimerized' in packing:
        p = 'dimerized'
        c = 'green'
    elif 'edge_brickwork' in packing:
        p = 'edge_brickwork'
        c = 'pink'
    elif 'edge_slipped_stack' in packing:
        p = 'edge_slipped_stack'
        c = 'cyan'
    elif 'edge_herringbone' in packing:
        p = 'edge_herringbone'
        c = 'yellow'
    elif 'sparse' in packing:
        p = 'sparse'
        c = 'gray'
    else:
        p = 'other'
        c = 'black'
    return p, c


empirical_list = ['3', '4', '5', '6', '7', '8', '9', '10', '11', '13', '14', '15', '18', '19', '20', '21', '22', '23', '24', '28', '29', '30', '31', '36', '37', '40', '44', '46', '48', '50', '53', '55', '56', '57', '61', '64', '65', '66', '67', '69', '71', '72', '74', '75', '77', '78', '80', '82', '83', '84', '86', '88', '89', '91', '93', '95', '96', '100', '101', '103', '108', '110', '116', '117', '119', '120', '128', '132', '133', '136', '144', '146', '147', '148', '156', '157', '161', '164', '166', '169', '170', '178', '180', '181', '182', '185', '188', '189', '194', '196', '197', '198', '199', '202', '203', '205']
cans = []
data = []
namelist = []
for fn in glob.glob('*.json'):
    d = parse_json(fn)
    packingdata = d['packing_data']
    name = d['cifname']
    p, c = emprical_packing(packingdata)

    # md = d['mol_data']
    # bonedata = md['bone']
    # bone_nrings = bonedata['nrings']
    # bone_length = bonedata['lp']
    # # bone_width = md['bone']['lq']
    # # bone_linearity = md['bone']['linearity']
    # # bone_can = md['bone']['canonical_smiles']
    # md['sidechains'] = sorted(md['sidechains'], key=lambda x:x['natoms'], reverse=True)[0]
    #
    # sc_branchrank = md['sidechains']['branch_rank']
    # sc_can = md['sidechains']['canonical_smiles']
    # sc_volume = md['sidechains']['volume']
    # sc_maxrank = md['sidechains']['maxrank']
    # sc_des_RadiusOfGyration = md['sidechains']['descriptors']['RadiusOfGyration']
    # sc_des_SpherocityIndex = md['sidechains']['descriptors']['SpherocityIndex']
    # sc_des_Asphericity = md['sidechains']['descriptors']['Asphericity']
    # sc_des_Eccentricity = md['sidechains']['descriptors']['Eccentricity']
    # sc_des_InertialShapeFactor = md['sidechains']['descriptors']['InertialShapeFactor']
    # sc_similarity_v_tips = md['sidechains']['similarity_v_tips']
    #
    # scud = md['sidechains']['scu']
    # scu_volume = scud['volume']
    # scu_can = scud['canonical_smiles']
    # scu_similarity_v_tipsu = scud['similarity_v_tipsu']
    # scu_des_RadiusOfGyration = scud['descriptors']['RadiusOfGyration']
    # scu_des_SpherocityIndex = scud['descriptors']['SpherocityIndex']
    # scu_des_Asphericity = scud['descriptors']['Asphericity']
    # scu_des_Eccentricity = scud['descriptors']['Eccentricity']
    # scu_des_InertialShapeFactor = scud['descriptors']['InertialShapeFactor']
    bulky_scs = 0
    for scd in d['mol_data']['sidechains']:
        if scd['natoms'] >= 3:
            bulky_scs += 1

    if bulky_scs == 2:
    # if bulky_scs == 2 and name in empirical_list:
        namelist.append(name)
        v = dict(
            similarity_v_tipsu=d['mol_data']['sidechains'][0]['scu']['similarity_v_tipsu'],
            volumeu=d['mol_data']['sidechains'][0]['scu']['volume'],
            Asphericityu=d['mol_data']['sidechains'][0]['scu']['descriptors']['Asphericity'],
            SpherocityIndexu=d['mol_data']['sidechains'][0]['scu']['descriptors']['SpherocityIndex'],
            RadiusOfGyrationu=d['mol_data']['sidechains'][0]['scu']['descriptors']['RadiusOfGyration'],
            InertialShapeFactoru=d['mol_data']['sidechains'][0]['scu']['descriptors']['InertialShapeFactor'],
            Eccentricityu=d['mol_data']['sidechains'][0]['scu']['descriptors']['Eccentricity'],

            similarity_v_tips=d['mol_data']['sidechains'][0]['similarity_v_tips'],
            volume=d['mol_data']['sidechains'][0]['volume'],
            Asphericity=d['mol_data']['sidechains'][0]['descriptors']['Asphericity'],
            SpherocityIndex=d['mol_data']['sidechains'][0]['descriptors']['SpherocityIndex'],
            RadiusOfGyration=d['mol_data']['sidechains'][0]['descriptors']['RadiusOfGyration'],
            InertialShapeFactor=d['mol_data']['sidechains'][0]['descriptors']['InertialShapeFactor'],
            Eccentricity=d['mol_data']['sidechains'][0]['descriptors']['Eccentricity'],

            color=c,
            packing=p,
        )
        sc_can = d['mol_data']['sidechains'][0]['canonical_smiles']
        # if d['mol_data']['bone']['nrings'] == 5:
        #     if sc_can not in cans:
        #         data.append(v)
        #         cans.append(sc_can)
        if d['mol_data']['bone']['nrings'] == 5:
            data.append(v)
            if p == 'brickwork':
                print(name, v['Eccentricityu'], p)

            # print(v)
# print(namelist)
print(len(cans))
print(len(set(cans)))
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.rc('font',family='Times New Roman')
plt.rcParams["font.size"] = 12

keys = list(data[0].keys())
keys = [ke for ke in keys if ke not in ['packing', 'color']]
possible_packings = list(set([d['packing'] for d in data]))


for i in range(len(keys)):
    k = keys[i]
    histd = []
    cols = []
    for packing_type in possible_packings:
        x = []
        for d in data:
            if d['packing'] == packing_type:
                x.append(d[k])
                col = d['color']

        histd.append(x)
        cols.append(col)
    bins = 20
    if k != 'SpherocityIndexu':
        continue
    plt.hist(histd, bins, color=cols, stacked=True, label=possible_packings)
    plt.legend()
    plt.tight_layout()
    plt.margins(0.05, 0.35)
    # plt.savefig('hist-{}.png'.format(k), dpi=600)
    plt.savefig('hist-{}.png'.format(k))
    plt.cla()


# for i in range(len(keys)):
#     k1 = 'packing'
#     k2 = keys[i]
#     px = [d[k2] for d in data]
#     py = [d[k1] for d in data]
#     c = [p2c(d['packing']) for d in data]
#     plt.scatter(px, py, c=c, alpha=0.2)
#     plt.savefig('{}-{}.png'.format(k1, k2))
#     plt.clf()

# for i in range(len(keys)):
#     for j in range(i, len(keys)):
#         k1 = keys[i]
#         k2 = keys[j]
#         px = [d[k1] for d in data]
#         py = [d[k2] for d in data]
#         c = [p2c(d['packing']) for d in data]
#         if k1=='packing':
#             plt.scatter(px, py, c=c, alpha=0.2)
#             plt.savefig('{}-{}.png'.format(k1, k2))
#             plt.clf()





# plt.scatter(d[:, 0].astype(np.float), d[:, 1].astype(np.float), c=d[0][3], alpha=float(d[0][4]), label=d[0][2])



