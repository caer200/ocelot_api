"""
test convex/concave hull method
"""
from ocelot.routines.geometry import alpha_shape, get_proj_point2plane, coord_transform
from ocelot.schema.omol import OMol
from descartes import PolygonPatch
import pylab as pl
from shapely.geometry import Polygon
import numpy as np

def plot_hull(xyzfile):
    omol = OMol.from_file(xyzfile)
    name = xyzfile[:-4]
    ref_proj_pts = [get_proj_point2plane(rs.coords, omol.backbone.vo_fit, omol.backbone.geoc) for rs in omol.backbone.msites]
    ref_2dpts = [coord_transform(omol.backbone.vp_fit, omol.backbone.vq_fit, omol.backbone.vo_fit, p3d)[:2] for p3d in ref_proj_pts]
    ref_2dpts = np.array(ref_2dpts)
    convex_hull = Polygon(ref_2dpts).convex_hull
    concave_hull, edges = alpha_shape(ref_2dpts, alpha=0.7)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    patch1 = PolygonPatch(concave_hull, fc='red', ec='#000000', fill=True, zorder=-1, alpha=0.2)
    patch2 = PolygonPatch(convex_hull, fc='blue', ec='#000000', fill=True, zorder=-1, alpha=0.2)
    ax.add_patch(patch1)
    ax.add_patch(patch2)
    ax.plot(ref_2dpts[:, 0], ref_2dpts[:, 1], 'ro')
    pl.savefig('hull_{}.png'.format(name))

plot_hull('2b.xyz')
plot_hull('3b.xyz')
plot_hull('6x.xyz')
plot_hull('adt.xyz')
