import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull


# 8 points defining the cube corners
pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0.5, 0.5, 0.5]])

hull = ConvexHull(pts)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Plot defining corner points
ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")

# 12 = 2 * 6 faces are the simplices (2 simplices per square face)
for s in hull.vertices:
    print(pts[s])
    # ax.plot(s, 'ro')
    # s = np.append(s, s[0])  # Here we cycle back to the first coordinate

# Make axis label
# for i in ["x", "y", "z"]:
#     [enter image description here][1]eval("ax.set_{:s}label('{:s}')".format(i, i))

plt.show()