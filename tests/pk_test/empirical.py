import matplotlib.pyplot as plt
import numpy as np

plt.switch_backend('agg')
plt.rc('font', family='Times New Roman')
plt.rcParams["font.size"] = 12

with open('packing.txt', 'r') as f:
    lines = f.readlines()

data = []
for l in lines[1:]:
    name, packingstring, lp, dsc = l.strip().split()
    lp = float(lp)
    dsc = float(dsc)
    packing = packingstring.split('-')
    if 'brickwork' in packing:
        p = 'brickwork'
        c = 'red'
        alpha = 1
    elif 'slipped_stack' in packing:
        p = 'slipped_stack'
        c = 'blue'
        alpha = 1
    elif 'herringbone' in packing:
        p = 'herringbone'
        # p = 'other'
        c = 'gold'
        alpha = 1
    elif 'dimerized' in packing:
        p = 'dimerized'
        # p = 'other'
        c = 'green'
        alpha = 1
    elif 'edge_brickwork' in packing:
        p = 'edge_brickwork'
        # p = 'brickwork'
        c = 'pink'
        alpha = 0.5
    elif 'edge_slipped_stack' in packing:
        p = 'edge_slipped_stack'
        # p = 'slipped_stack'
        c = 'cyan'
        alpha = 0.5
    elif 'edge_herringbone' in packing:
        p = 'edge_herringbone'
        # p = 'herringbone'
        # p = 'other'
        c = 'yellow'
        alpha = 0.5
    elif 'sparse' in packing:
        p = 'sparse'
        # p = 'other'
        c = 'gray'
        alpha = 0.5
    else:
        p = 'other'
        c = 'black'
        alpha = 0.2
    # data.append([lp, dsc, p, c, alpha])
    data.append([lp / (dsc - 1.8), p, c, name])

empirical_list = []
for i in data:
    if 2.4 > i[0] > 1.7:
        empirical_list.append(i[3])
print(empirical_list)


possible_packings = set([d[1] for d in data])
possible_packings = list(possible_packings)

histd = []
colors = []
labels = []
twobulkysc = ['3', '118', '40', '205', '177', '11', '151', '49', '140', '120', '136', '184', '82', '156', '114', '68',
              '174', '33', '89', '16', '104', '58', '19', '123', '99', '71', '188', '197', '39', '128', '46', '76',
              '54', '107', '190', '101', '164', '43', '63', '86', '173', '57', '186', '67', '113', '147', '133', '176',
              '55', '121', '117', '207', '167', '199', '37', '148', '30', '83', '196', '27', '100', '112', '152', '181',
              '171', '134', '59', '162', '119', '204', '132', '183', '77', '28', '66', '115', '202', '32', '135', '103',
              '93', '60', '200', '47', '91', '21', '175', '179', '161', '129', '111', '194', '4', '45', '6', '125',
              '131', '53', '2', '180', '88', '153', '158', '139', '62', '41', '126', '61', '75', '155', '185', '38',
              '52', '182', '50', '124', '149', '85', '1', '122', '73', '159', '84', '51', '157', '144', '36', '80',
              '105', '44', '178', '64', '102', '193', '187', '191', '18', '146', '17', '87', '165', '74', '138', '108',
              '69', '154', '48', '143', '116', '10', '34', '70', '198', '195', '169', '110', '5', '150', '201', '35']
for packing in possible_packings:
    labels.append(packing)
    thishist = []
    cc = ""
    for d in data:
        if d[1] == packing and d[3] in twobulkysc:
            thishist.append(d[0])
            cc = d[2]
    histd.append(thishist)
    colors.append(cc)

bins = 40
plt.hist(histd, bins, color=colors, stacked=True, label=labels)
plt.ylabel('Counts')
plt.xlabel('L/D')
plt.legend()
plt.margins(0.05)
plt.tight_layout()
plt.savefig("hist.png", dpi=600)
