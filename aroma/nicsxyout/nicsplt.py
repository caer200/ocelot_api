import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend('agg')
plt.rc('font',family='Times New Roman')
plt.rcParams["font.size"] = 14

spicies = ['Br', 'Cl', 'CN', 'CNBr', 'F', 'H']
states = ['s0', 't1_delta']

sigmafiles = glob.glob('*sigma.zz')
datafiles = glob.glob('*.data')
totalfiles = glob.glob('*total.zz')

cc = ['r', 'b', 'k']

for sp in spicies:
    for state in states:
        for sf in sigmafiles:
            if sp in sf and state in sf:
                sigma = sf
                break
        for sf in totalfiles:
            if sp in sf and state in sf:
                total = sf
                break
        for sf in datafiles:
            if sp in sf and state in sf:
                data = sf
                break

        sigmavalues = np.loadtxt(sigma)[:, 1]
        totalvalues = np.loadtxt(total)[:, 1]
        d = np.loadtxt(data, skiprows=1)
        x = d[:, 0]
        with open(data, 'r') as f:
            xticks = f.readlines()[0].strip().split()
        xticks = [float(t) for t in xticks]
        plt.scatter(x, totalvalues-sigmavalues, c='b')
        plt.scatter(x, totalvalues, c='k')
        plt.scatter(x, sigmavalues, c='r')
        vc = ['gray', 'black']
        i=1
        for xv in xticks:
            plt.axvline(xv, c=vc[i])
            i = 1-i
        plt.xlim([x[0], x[-1]])
        plt.ylim([-22, 5])
        plt.xlabel(r'scan path ($\rm{\AA}$)')
        plt.ylabel('NICSzz (ppm)')
        plt.savefig(sp+'-'+state + '.png', dpi=800)
        plt.tight_layout()
        plt.clf()


