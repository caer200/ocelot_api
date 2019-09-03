import numpy as np
import sys
import math


class LoopSearcher:

    def __init__(self, nbmap):
        """
        smallest set of smallest rings
        giving a connection table, find all possible loops with a certain loop size
        nbmap[i] does not contain i itself
        check 10.1073pnas.0813040106 for a better solution
        :param nbmap: connection table
        """
        self.nbmap = nbmap

    def expand(self, path_set):
        """
        the path will never intersect itself
        """
        new_path_set = []

        if len(path_set) == 1 and len(path_set[0]) == 1:
            start = path_set[0][0]
            for n in self.nbmap[start]:
                if n != start:
                    new_path_set.append([start, n])
            return new_path_set

        for path in path_set:
            for n in self.nbmap[path[-1]]:
                if n not in path:
                    new_path_set.append(path + [n])

        if len(new_path_set) == 0:
            return path_set
        else:
            return new_path_set

    def alex_method(self, loopsize):
        """
        I figured this out but I'm sure I'm not the first
        :param loopsize:
        :return:
        """
        all_idx = list(range(len(self.nbmap)))
        loop_found = []
        visited = []
        while len(visited) != len(all_idx):
            unvisited = [i for i in all_idx if i not in visited]
            start = unvisited[0]
            path_set = [[start]]
            for i in range(loopsize - 1):
                path_set = self.expand(path_set)
            for path in path_set:
                if path[0] in self.nbmap[path[-1]]:
                    if set(path) not in [set(loop) for loop in loop_found]:
                        loop_found.append(path)
                        visited += [p for p in path if p not in visited]
                    break
            if start not in visited:
                visited += [start]
        return loop_found


class Geo:

    @staticmethod
    def angle_btw(v1, v2, output='radian'):
        """
        get angle between two vectors
        :param v1:
        :param v2:
        :param output:
        :return:
        """
        v1_u = Geo.unify(v1)
        v2_u = Geo.unify(v2)
        angle_radian = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        if output == 'radian':
            return angle_radian
        elif output == 'degree':
            return np.rad2deg(angle_radian)
        else:
            sys.exit('output type not recongized')

    @staticmethod
    def norm(v):
        l = 0
        if isinstance(v, float) or isinstance(v, int):
            return v
        for i in range(len(v)):
            l += v[i]**2
        return math.sqrt(l)

    @staticmethod
    def unify(v):
        return v/Geo.norm(v)

    @staticmethod
    def dist_pt2line(p, q, r):
        """
        https://stackoverflow.com/questions/50727961/
        :param p: end of segment
        :param q: end of segment
        :param r: the point
        :return:
        """
        x = p-q
        t = np.dot(r-q, x)/np.dot(x, x)
        return np.linalg.norm(t*(p-q)+q-r)

    @staticmethod
    def genlinepts(a, b, stepsize):
        pts = []
        for i in range(int(np.floor(np.linalg.norm(b-a)/stepsize))+1):
            pts.append(a + stepsize*Geo.unify(b-a)*i)
        return pts

class GeoFitter3d:
    """
    https://www.ltu.se/cms_fs/1.51590!/svd-fitting.pdf
    https://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
    https://math.stackexchange.com/questions/2378198/computing-least-squares-error-from-plane-fitting-svd
    https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """
    def __init__(self, pts):
        self.pts = np.zeros((len(pts), 3))
        for i in range(len(pts)):
            for j in range(3):
                self.pts[i][j] = pts[i][j]

    @staticmethod
    def iscollinear(pts, tol=1e-3):
        v, m, e = GeoFitter3d.linear_fit(pts)
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
        vector = Geo.unify(vt[0])
        error = s[1]**2 + s[2]**2
        return vector, ptsmean, error

    @staticmethod
    def plane_fit(pts):
        if len(pts) < 3:
            sys.exit('less than 3 pts to fit a plane!')
        if GeoFitter3d.iscollinear(pts):
            sys.exit('trying to fit a plane with 3 collinear points!')
        ptsmean = np.mean(pts, axis=0)
        u, s, vt = np.linalg.svd(pts - ptsmean)
        normal = Geo.unify(vt[2])
        error = s[2]**2
        return normal, ptsmean, error

