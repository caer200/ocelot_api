from pymatgen.io.gaussian import GaussianOutput
from pymatgen.util.plotting import pretty_plot
from matplotlib.mlab import normpdf
import numpy as np
import scipy.constants as cst

def pltuvvis(fn='tdsinglet.log'):
    go = GaussianOutput(fn)
    plt = pretty_plot(12, 8)

    sigma = 0.05
    step = 0.01
    transitions = go.read_excitation_energies()

    minval = min([val[0] for val in transitions]) - 5.0 * sigma
    maxval = max([val[0] for val in transitions]) + 5.0 * sigma
    npts = int((maxval - minval) / step) + 1

    eneval = np.linspace(minval, maxval, npts)  # in eV
    lambdaval = [cst.h * cst.c / (val * cst.e) * 1.e9
                 for val in eneval]  # in nm

    # sum of gaussian functions
    spectre = np.zeros(npts)
    for trans in transitions:
        spectre += trans[2] * normpdf(eneval, trans[0], sigma)
    spectre /= spectre.max()
    plt.plot(lambdaval, spectre, "r-", label="Simulated Spectrum")

    data = {"energies": eneval, "lambda": lambdaval, "xas": spectre}
    plt.xlim([min(lambdaval), max(lambdaval)])
    # plot transitions as vlines
    plt.vlines([val[1] for val in transitions],
               0.,
               [val[2] for val in transitions],
               color="blue",
               label="Eigenvalues",
               linewidth=2)

    plt.xlabel("$\\lambda$ (nm)")
    plt.ylabel("$f$ (a.u.)")
    plt.legend(prop={'size': 22})
    plt.tight_layout()
    plt.savefig('uvvis.png')