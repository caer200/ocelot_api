
def parse(logname):
    with open(logname, 'r') as f:
        lines = f.readlines()

    ns = 0
    for i in range(len(lines)):
        if 'Charge' in lines[i] and 'Multiplicity' in lines[i]:
            cursor = i+1
            while len(lines[cursor].strip().split()) == 4:
                ns += 1
                cursor += 1
            break

    tensor_zzs = []
    for i in range(len(lines)):
        if 'GIAO Magnetic shielding tensor' in lines[i]:
            for j in range(i+1, i+1+ns*5):
                if 'Bq' in lines[j]:
                    tensor_zzs.append(-float(lines[j+3].strip().split()[-1]))
            break
    with open(logname[:-4]+'.zz', 'w') as f:
        for i in range(len(tensor_zzs)):
            f.write(str(i) + '\t' + str(tensor_zzs[i]) + '\r\n')


import glob
for i in glob.glob('*.log'):
    parse(i)