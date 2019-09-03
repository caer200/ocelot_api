import math
import sys
import warnings
# import sympy
# import sympy.geometry

import numpy as np


def norm(v):
    l = 0
    if isinstance(v, float) or isinstance(v, int):
        return v
    for i in range(len(v)):
        l += v[i] ** 2
    return math.sqrt(l)


def are_collinear(v1, v2):
    return abs((np.array(v1) @ np.array(v2) / (norm(v1) * norm(v2))) - 1) < 1e-5


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
        sys.exit('output type not recongized')


def unify(v):
    return v / norm(v)


def dist_pt2line(p, q, r):
    """
    https://stackoverflow.com/questions/50727961/
    :param p: end of segment
    :param q: end of segment
    :param r: the point
    :return:
    """
    x = p - q
    t = np.dot(r - q, x) / np.dot(x, x)
    return np.linalg.norm(t * (p - q) + q - r)


def genlinepts(a, b, stepsize):
    pts = []
    for i in range(int(np.floor(np.linalg.norm(b - a) / stepsize)) + 1):
        pts.append(a + stepsize * unify(b - a) * i)
    return pts


def coord_transform(p, q, o, coords):
    """
    get the new coord under kart system defined by p, q, o as mutual ort unit vectors
    :param p:
    :param q:
    :param o:
    :param coords:
    :return:
    """
    newxyz = np.array([p, q, o]).T
    return np.linalg.solve(newxyz, coords)


def coord_reverse(p, q, o, coords):
    mat = np.array([p, q, o])
    return np.matmul(coords, mat.T)


def cart2frac(coords, pbc):
    """
    convert to frac
    :param coords:
    :param pbc: a 3x3 mat, [a, b, c]
    :return:
    """
    return np.matmul(coords, np.linalg.inv(pbc))


def frac2cart(coords, pbc):
    return np.matmul(coords, pbc)


def cart2polar(coords):
    ptsnew = np.zeros(coords.shape)
    xy = coords[0] ** 2 + coords[1] ** 2
    ptsnew[0] = np.sqrt(xy + coords[2] ** 2)
    ptsnew[1] = np.arctan2(np.sqrt(xy), coords[2])  # for elevation angle defined from Z-axis down
    ptsnew[2] = np.arctan2(coords[1], coords[0])
    return ptsnew


def abcabg2pbc(abcabg):
    a, b, c, alpha, beta, gamma = abcabg
    cosdelta = (math.cos(math.radians(alpha)) - math.cos(math.radians(beta)) * math.cos(math.radians(gamma))) / \
               (math.sin(math.radians(beta)) * math.sin(math.radians(gamma)))
    sindelta = math.sqrt(1 - cosdelta ** 2)
    pbc = np.zeros((3, 3))
    pbc[0] = a * np.array([1.0, 0.0, 0.0])
    pbc[1] = b * np.array([math.cos(math.radians(gamma)), math.sin(math.radians(gamma)), 0.0])
    pbc[2] = c * np.array([math.cos(math.radians(beta)),
                           math.sin(math.radians(beta)) * cosdelta, math.sin(math.radians(beta)) * sindelta])
    return pbc

# def get_sympypt_coords(pt, t='pi/2'):
#     """
#     somehow the point you get from plane.arbitrary_point can not be directly evaluated
#     e.g.
#
#     import sympy.geometry as geo
#     a = geo.Point(1, 2, 3)
#     b = geo.Point(1, 1, 2)
#     c = geo.Point(1, 22, 32)
#     p = geo.Plane([a, b, c])
#     pt = plane.arbitrary_point
#     pt.evalf(subs={'t': 2})  # does not work
#     pt[1].evalf(subs={'t': 2})  # does not work
#
#     it looks like the 't' in the expression is not recognized, this is just silly af
#     :param pt:
#     :return:
#     """
#     return [sympy.S(str(c)).evalf(subs={'t': t}) for c in pt]


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
    if axis is 'x':
        axis = np.array([1, 0, 0])
    elif axis is 'y':
        axis = np.array([0, 1, 0])
    elif axis is 'z':
        axis = np.array([0, 0, 1])
    else:
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def get_plane_param(normal, pt):
    """
    ax + by + cz + d = 0
    :param normal:
    :param pt:
    :return:
    """
    a = normal[0]
    b = normal[1]
    c = normal[2]
    d = (-1) * (a * pt[0] + b * pt[1] + c * pt[2])
    return np.array([a, b, c, d])

def arbitrary_noraml(v):
    if v[1] == 0.0 and v[2] == 0:
        if v[0] == 0.0:
            raise ValueError('zero vector')

        return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])


def get_proj_point2plane(pt, normal, ptonplane):
    """
    :param pt:
    :param normal:
    :param ptonplane:
    :return:
    """
    v_ptonplane2pt = pt - ptonplane
    angle1 = angle_btw(v_ptonplane2pt, normal)
    if angle1 == 0:
        return ptonplane
    aus = norm(v_ptonplane2pt) * math.cos(angle1) * normal
    v_ptonplane2pedal = v_ptonplane2pt - aus
    return ptonplane + v_ptonplane2pedal


def alpha_shape(points, alpha=0.7):
    """
    0.7 seems work for adt
    https://gist.github.com/dwyerk/10561690
    Compute the alpha shape (concave hull) of a set
    of points.
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    from shapely.ops import cascaded_union, polygonize
    from scipy.spatial import Delaunay
    import shapely.geometry as geometry
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull

    # coords = np.array([point[0] for point in points])
    coords = np.array(points)
    tri = Delaunay(coords)
    triangles = coords[tri.simplices]
    a = ((triangles[:,0,0] - triangles[:,1,0]) ** 2 + (triangles[:,0,1] - triangles[:,1,1]) ** 2) ** 0.5
    b = ((triangles[:,1,0] - triangles[:,2,0]) ** 2 + (triangles[:,1,1] - triangles[:,2,1]) ** 2) ** 0.5
    c = ((triangles[:,2,0] - triangles[:,0,0]) ** 2 + (triangles[:,2,1] - triangles[:,0,1]) ** 2) ** 0.5
    s = ( a + b + c ) / 2.0
    areas = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    circums = a * b * c / (4.0 * areas)
    filtered = triangles[circums < (1.0 / alpha)]
    edge1 = filtered[:,(0,1)]
    edge2 = filtered[:,(1,2)]
    edge3 = filtered[:,(2,0)]
    edge_points = np.unique(np.concatenate((edge1,edge2,edge3)), axis = 0).tolist()
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

class Fitter:
    """
    default for 3d
    https://www.ltu.se/cms_fs/1.51590!/svd-fitting.pdf
    https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    https://math.stackexchange.com/questions/2378198/computing-least-squares-error-from-plane-fitting-svd
    https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """

    # def __init__(self, pts):
    #     self.pts = np.zeros((len(pts), len(pts[0])))
    #     for i in range(len(pts)):
    #         for j in range(len(pts[0])):
    #             self.pts[i][j] = pts[i][j]

    @staticmethod
    def iscollinear(pts, tol=1e-3):
        v, m, e = Fitter.linear_fit(pts)
        if e > tol:
            return 0
        return 1

    @staticmethod
    def linear_fit(pts):
        """
        vector * t + ptsmean = pt on fitted line
        :param pts: a (n, 3) array
        :return:
        """
        if len(pts) < 2:
            sys.exit('less than 2 pts to fit a line!')
        ptsmean = np.mean(pts, axis=0)
        u, s, vt = np.linalg.svd(pts - ptsmean)
        vector = unify(vt[0])
        error = s[1] ** 2 + s[2] ** 2
        return vector, ptsmean, error

    @staticmethod
    def plane_fit(pts):
        """
        this only returns a normal, we do not have a plane equation here
        :param pts:
        :return:
        """
        if len(pts) < 3:
            sys.exit('less than 3 pts to fit a plane!')
        if Fitter.iscollinear(pts):
            warnings.warn('trying to fit a plane with collinear points!')
        ptsmean = np.mean(pts, axis=0)
        u, s, vt = np.linalg.svd(pts - ptsmean)
        normal = unify(vt[2])
        error = s[2] ** 2
        return normal, ptsmean, error
