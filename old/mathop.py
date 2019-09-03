import math
import sys
import itertools
import warnings
import numpy as np
import operator
from scipy.optimize import leastsq


def plane_fitting(pts):
    """
    TODO profiling
    :param pts: (n, 3) np array
    :return:
    """
    if len(pts) < 3:
        warnings.warn('less than 3 pts to fit a plane!')
    elif len(pts) == 3:
        if iscolinear(pts[0], pts[1], pts[2]):
            warnings.warn('plane project failed as 3 pts on a line')
            sys.exit(1)
        else:
            param, norm = get_plane_param(pts[0], pts[1], pts[2])
            return param, norm
    else:
        # TODO you need another if to rule out collinear situation for many points
        # https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
        init_pts = np.zeros((3, 3), dtype=float)
        for idx_comb in itertools.combinations(range(len(pts)), 3):
            if not iscolinear(pts[idx_comb[0]], pts[idx_comb[1]], pts[idx_comb[2]]):
                init_pts[0] = pts[idx_comb[0]]
                init_pts[1] = pts[idx_comb[1]]
                init_pts[2] = pts[idx_comb[2]]
                break
        p0, init_norm = get_plane_param(init_pts[0], init_pts[1], init_pts[2])

        def f_min(xmat, p):
            plane_xyz = p[0:3]
            distance = (plane_xyz * xmat.T).sum(axis=1) + p[3]
            return distance / np.linalg.norm(plane_xyz)

        def residuals(paramas, signal, xmat):
            return f_min(xmat, paramas)

        sol = leastsq(residuals, np.array(p0), args=(None, pts.T))
        fin = np.array(sol[0] / sol[0][3])
        return fin, unify(fin[:3])


def min_widx(values):
    min_index, min_value = min(enumerate(values), key=operator.itemgetter(1))
    return min_index, min_value


def max_widx(values):
    max_index, max_value = max(enumerate(values), key=operator.itemgetter(1))
    return max_index, max_value


def max_distance(pts):
    pts = np.array(pts)
    distances = []
    for pair_idx in itertools.combinations(range(len(pts)), 2):
        i, j = pair_idx
        v = pts[i] - pts[j]
        distances.append((np.linalg.norm(v), v))
    distances = sorted(distances, key=lambda x: x[0], reverse=True)
    return distances[0]


def linear_fit(xys):
    """
    ax+b = y
    :param xys:
    :return:
    """
    x = xys[:, 0]
    y = xys[:, 1]
    out = np.polyfit(x, y, 1, full=True)
    a = out[0][1]
    b = out[0][0]
    residual = out[1][0]
    return a, b, residual


def loop_finding(nblist, loopsize=6):
    """
    giving a nrnbr table, find all possible loops with a certain loop size
    nrnbr[i] does not contain i
    :param nblist:
    :param loopsize: size of the ring
    :return:
    """
    all_idx = list(range(len(nblist)))
    loop_found = []
    visited = []
    while len(visited) != len(all_idx):
        unvisited = [i for i in all_idx if i not in visited]
        start = unvisited[0]
        path_set = [[start]]
        for i in range(loopsize - 1):
            path_set = expand_path_set(nblist, path_set)
        for path in path_set:
            if path[0] in nblist[path[-1]]:
                if set(path) not in [set(loop) for loop in loop_found]:
                    loop_found.append(path)
                    visited += [p for p in path if p not in visited]
                break
        if start not in visited:
            visited += [start]
    return loop_found

def linear_fit_3d(pts):
    """
    https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    https://math.stackexchange.com/questions/2378198/computing-least-squares-error-from-plane-fitting-svd
    vector * t + ptsmean = pt on fitted line
    :param pts: a (n, 3) array
    :return:
    """
    ptsmean = pts.mean(axis=0)
    uu, dd, vv = np.linalg.svd(pts - ptsmean)
    vector = unify(vv[0])
    error = abs(dd[-1])
    return vector, ptsmean, error

def expand_path_set(nblist, path_set):
    """
    the path will never intersect itself
    :param nblist:
    :param path_set:
    :return:
    """
    new_path_set = []

    if len(path_set) == 1 and len(path_set[0]) == 1:
        start = path_set[0][0]
        for next in nblist[start]:
            if next != start:
                new_path_set.append([start, next])
        return new_path_set

    for path in path_set:
        for next in nblist[path[-1]]:
            if next not in path:
                new_path_set.append(path + [next])

    if len(new_path_set) == 0:
        return path_set
    else:
        return new_path_set


def iscolinear(p1, p2, p3):
    # TODO add support for many points
    v1 = np.array(p1) - np.array(p3)
    v2 = np.array(p1) - np.array(p2)
    if angle_btw(v1, v2) < 0.00001:
        return True
    else:
        return False


def rotate_along_2pts(s1, s2, theta, s0, thetaunit='degree'):
    """
    this can actually be achieved with rotate_along_axis
    but I coded this earlier so I still keep it

    ax2 s2
    |
    |
    |3--------d-------site_to_be_rotated s0         y
    |                                               |
    ax1 s1                                          z----x

    :param s1:
    :param s2:
    :param theta:
    :param s0:
    :param thetaunit:
    """
    if thetaunit == 'degree':
        theta = np.deg2rad(theta)

    if np.allclose(s0, s1) or np.allclose(s0, s2):
        return np.array(s0)

    else:
        # trans
        s0_trans = kart_trans(s1, s2, s0, s0)
        s1_trans = kart_trans(s1, s2, s0, s0)
        s2_trans = kart_trans(s1, s2, s0, s0)
        # find foot point
        s3_coords = s1_trans + np.dot((s0_trans - s1_trans),
                                      (s2_trans - s1_trans)) * \
            unify(s2_trans - s1_trans) / np.linalg.norm((s2_trans - s1_trans))
        d = np.linalg.norm(s3_coords - s0_trans)
        # initialize theta
        s0_theta = np.arctan((s0_trans[2] - s3_coords[2]) / (s0_trans[0] - s3_coords[0]))

        if abs(s0_theta) > 0.00001:
            warnings.warn('theta initialization failed')
            sys.exit(1)

        if (s0_trans[0] - s3_coords[0]) > 0:
            s0_new_theta = s0_theta + theta
            s0_new_coords_x = d * np.cos(s0_new_theta) + s3_coords[0]
        else:
            s0_new_theta = s0_theta - theta
            s0_new_coords_x = -d * np.cos(s0_new_theta) + s3_coords[0]

        s0_new_coords_y = s0_trans[1]
        s0_new_coords_z = d * np.sin(s0_new_theta) + s3_coords[2]
        s0_new_coords = np.array([s0_new_coords_x, s0_new_coords_y, s0_new_coords_z])
        return s0_new_coords


def kart_reverse(p1, p2, p3, newcoord):
    """
    be aware that the ps here should be identical to the ones used in trans
    :param p1:
    :param p2:
    :param p3:
    :param newcoord:
    :return:
    """
    trans_matrix = plane_kart_3pts(p1, p2, p3).T
    mat = np.matrix(trans_matrix)
    product = np.matmul(mat, newcoord)
    oldcoord = np.squeeze(np.asarray(product))
    return oldcoord


def kart_trans(p1, p2, p3, coord):
    """
    get the new coord under kart system defined by the plane of p1, p2, p3
    p1, p2, p3 are defined in real kart
    :param p1:
    :param p2:
    :param p3:
    :param coord:
    :return:
    """
    trans_matrix = plane_kart_3pts(p1, p2, p3).T
    old_coords = np.array(coord)
    return np.linalg.solve(trans_matrix, old_coords)


def plane_kart_3pts(p1, p2, p3):
    """
    2     ^
      3   y  x>   z.
    1
    :param p1:
    :param p2:
    :param p3:
    :return:
    """
    coord1 = np.array(p1)
    coord2 = np.array(p2)
    coord3 = np.array(p3)
    param, norm = get_plane_param(coord1, coord2, coord3)
    y_u = unify(coord2 - coord1)
    z_u = unify(norm)
    x_u = unify(np.cross(y_u, z_u))
    return np.array([x_u, y_u, z_u])


def get_plane_param(coord1, coord2, coord3):
    """
    ax + by + cz + d = 0
    :param coord1:
    :param coord2:
    :param coord3:
    :return:
    """
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)
    coord3 = np.array(coord3)
    v1 = coord1 - coord2
    v2 = coord3 - coord2
    normal = np.cross(v1, v2)
    a = normal[0]
    b = normal[1]
    c = normal[2]
    d = (-1) * (a * coord1[0] + b * coord1[1] + c * coord1[2])
    norm = np.array([a, b, c])
    return [a, b, c, d], norm


def rotate_along_axis(v, axis, theta, thetaunit='degree'):
    """
    rotate v along axis counterclockwise
    :param v:
    :param axis:
    :param theta:
    :param thetaunit:
    :return:
    """
    return np.dot(rotation_matrix(axis, theta, thetaunit), v)


def rotation_matrix(axis, theta, thetaunit='degree'):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta.
    np.dot(rotation_matrix(axis,theta_d), v)
    https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    :param axis: a list of 3 floats
    :param theta
    :param thetaunit
    """
    if thetaunit == 'degree':
        theta = np.deg2rad(theta)
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def fd_reci_2ndder(y, x, x0_index, step=1, accuracy='high'):
    """
    https://en.wikipedia.org/wiki/Finite_difference_coefficient
    x should be uniformly sampled over a seg
    :param y:
    :param x: length should be = (order+2)
    :param x0_index: index of x0 at which derivative is calculated
    :param step: step size based on index
    :param accuracy:
    :return: 1/y''
    """
    if accuracy == 'low':
        cfb = [1.0, -2.0, 1.0]
        cc = [1.0, -2.0, 1.0]
    else:
        cfb = [2.0, -5.0, 4.0, -1.0]  # forward and backward
        # cfb = [35.0/12.0, -26.0/3.0, 19.0/2.0, -14.0/3.0, 11.0/12.0]
        cc = [-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0]
    if len(x) < 2 * len(cfb):
        warnings.warn('too few points for finite diff at accuracy = ' + accuracy)
        return 10000.0
    summ = 0.0
    if x0_index == 0 or x0_index == 1 or x0_index == 2:  # forward fd
        for k in range(len(cfb)):
            summ += cfb[k] * y[k * step]
    elif x0_index == len(x) - 1 or x0_index == len(x) - 2 or x0_index == len(x) - 3:
        for k in range(len(cfb)):
            summ += cfb[k] * y[-1 + (-1 * k * step)]
    elif (x0_index + 2 * step) in range(len(x)) and (x0_index - 2 * step) in range(len(x)):
        for k in range(len(cc)):
            summ += cc[k] * y[int(round(x0_index + k * step - ((len(cc) - 1) / 2) * step))]
    else:
        warnings.warn('center fd failed as this kpt is too close to seg boundary')
        return 100000.0
    der2 = summ / ((x[0] - x[step]) ** 2)
    return 1.0 / der2


def ra_to_rb(ra_list):
    rb_list = []
    for i in ra_list:
        rb = float(i) * np.pi * 2 * (1.0 / 1.8897259885789)  # so vasp outcar kpt cart is 2pi/A ???
        rb_list.append(rb)
    return rb_list


def ev_to_ha(eigen_list):
    ha_list = []
    for i in eigen_list:
        ha = float(i) * 0.0367493
        ha_list.append(ha)
    return ha_list


def coords_converter2cart(original_coords, current_basis):
    """
    convert to cartesian coords
    :param original_coords: coords in certain basis
    :param current_basis: basis vectors in cartesian
    :return: [x, y, z]
    """
    coords = np.array(original_coords, dtype='float')
    basis = np.array(current_basis, dtype='float')
    r = np.matmul(coords, basis)
    return r


def unify(vector):
    """
    normalize vector
    :param vector:
    :return:
    """
    if np.linalg.norm(vector) == 0.0:
        return np.zeros((len(vector)))
    else:
        return vector / np.linalg.norm(vector)


def angle_btw(v1, v2, output='radian'):
    """
    get angle between two vectors
    :param v1:
    :param v2:
    :param output:
    :return:
    """
    v1_u = unify(v1)
    v2_u = unify(v2)
    angle_radian = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if output == 'radian':
        return angle_radian
    elif output == 'degree':
        return np.rad2deg(angle_radian)
    else:
        warnings.warn('output type not recongized')
        sys.exit(1)


def latparam2matrix(latparam):
    """
    [a, b, c, alpha, beta, gamma] --> [va, vb, vc]
    :param latparam:
    :return:
    """
    a, b, c, alpha, beta, gamma = [float(i) for i in latparam]
    cosdelta = (math.cos(math.radians(alpha)) - math.cos(math.radians(beta)) * math.cos(math.radians(gamma))) / \
               (math.sin(math.radians(beta)) * math.sin(math.radians(gamma)))
    sindelta = math.sqrt(1 - cosdelta ** 2)
    va = a * np.array([1.0, 0.0, 0.0])
    vb = b * np.array([math.cos(math.radians(gamma)), math.sin(math.radians(gamma)), 0.0])
    vc = c * np.array([math.cos(math.radians(beta)), math.sin(math.radians(beta)) * cosdelta,
                       math.sin(math.radians(beta)) * sindelta])
    volume = np.dot(va, np.cross(vb, vc))
    return np.array([va, vb, vc]), volume
