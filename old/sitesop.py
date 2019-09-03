import math
import itertools
import constants
import numpy as np
from scipy.spatial.distance import pdist, squareform

'''
operations for a list of msites
'''


def nrnbrmap(sites):
    bmat = bondmat(sites)
    ma = {}
    numsites = len(bmat)
    for i in range(numsites):
        nbs = []
        for j in range(numsites):
            if bmat[i][j] and j != i:
                nbs.append(j)  # TODO profile
        ma[i] = nbs
    return ma


def coordmat(sites):
    """
    coordinates matrix
    :param sites
    :return:
    """
    mat = np.empty((len(sites), 3))
    for i in range(len(sites)):
        for j in range(3):
            mat[i][j] = sites[i].coords[j]
    return mat


def bondmat(sites, distmat=None, co=1.3):
    """
    Bij = whether there is a bond between si and sj
    i is NOT bonded with itself
    :param sites
    :param distmat
    :param co: coefficient for cutoff
    :return:
    """
    if distmat is None:
        comat = coordmat(sites)
        distmat = squareform(pdist(comat))
    numsites = len(sites)
    mat = np.ones((numsites, numsites), dtype=bool)
    for i in range(numsites):
        mat[i][i] = False
        for j in range(i + 1, numsites):
            if hasattr(sites[i], 'symbol'):
                if distmat[i][j] > (constants.radius[sites[i].symbol] + constants.radius[sites[j].symbol]) * co:
                    mat[i][j] = False
                    mat[j][i] = False
            elif hasattr(sites[i], 'species_string'):
                if distmat[i][j] > (
                        constants.radius[sites[i].species_string] + constants.radius[sites[j].species_string]) * co:
                    mat[i][j] = False
                    mat[j][i] = False
    return mat


def sites_toxyz(sites, xyzname):
    with open(xyzname, 'w') as f:
        f.write(str(len(sites)) + '\r\n\r\n')
        for s in sites:
            f.write(s.symbol + ' ' + ' '.join([str(c) for c in s.coords]) + '\r\n')


def volume(sites):
    """
    # TODO profile, linalg inefficient
    http://wiki.bkslab.org/index.php/Calculate_volume_of_the_binding_site_and_molecules
    First, Lay a grid over the spheres.
    Count the number or points contained in the spheres (Ns).
    Count the number of points in the grid box (Ng).
    Calculate the volume of the grid box (Vb).
    :return:
    """
    box_density = 0.5  # \AA^3
    mat = np.empty((len(sites), 4))
    for i in range(len(sites)):
        for j in range(3):
            mat[i][j] = sites[i].coords[j]
        mat[i][3] = constants.van_radius[sites[i].symbol]

    box_min = math.floor(min(itertools.chain.from_iterable(mat))) - 2
    box_max = math.ceil(max(itertools.chain.from_iterable(mat))) + 2
    axis = np.arange(box_min, box_max + box_density, box_density)
    grid = np.array(np.meshgrid(axis, axis, axis)).T.reshape(-1, 3)
    ngps_total = len(grid)
    ngps_occu = 0
    for igp in range(len(grid)):
        for iap in range(len(mat)):
            dist = np.linalg.norm(grid[igp] - mat[iap][:3])
            if dist < mat[iap][3]:
                ngps_occu += 1
                break
    v = (ngps_occu / ngps_total) * ((box_max - box_min) ** 3)
    return v


def nelectrons(sites):
    ne = 0
    for s in sites:
        ne += s.properties['number']
    return ne
