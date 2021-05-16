import numpy as np
import ast
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay
import sys


v = (0.35, 0.25)
a = (0.4, 0.3)
b = (0.35, 0.35)
center = (0.45, 0.35)
points = [ (0.35, 0.35), (0.55, 0.15), (0.55, 0.25), ]
voronoi = [ (0.4, 0.2), (0.4, 0.2), (0.4, 0.2), ]
radii = [ 0.158114, 0.158114, 0.158114, ]


points = np.array(points)
voronoi = np.array(voronoi)
center = np.array(center)
v = np.array(v)

# a = (center+v)/2
# b = a+np.array([(v-center)[1], -(v-center)[0]])

fig = plt.figure()
ax = fig.add_subplot(111)

# tmp = (v-center)*100
# ax.plot([0,tmp[0]], [0,tmp[1]])
# ax.plot([0,tmp[1]], [0,-tmp[0]])
# plt.scatter(0, 0)


# plt.show()

# sys.exit()
if (a[0]-b[0]!=0):
    h = lambda x: (a[1]-b[1])/(a[0]-b[0])*(x-b[0])+b[1]
    ax.plot(np.linspace(0,0.5), h(np.linspace(0,0.5)))



f = lambda x: (center[1]-v[1])/(center[0]-v[0])*(x-v[0])+v[1]
#(-2.425, 0.35) 
# print(h(0.03)<=0.35)



# ax.plot(np.linspace(0,0.5), f(np.linspace(0,0.5)))

ax.plot(points[:,0], points[:,1], '-o', label="pnts", color="lime")
ax.plot(voronoi[:,0], voronoi[:,1], '-o', label="voronoi", color="orange")
ax.scatter(center[0], center[1], label="center", color="red")
ax.scatter(v[0], v[1], label="candidate", color="blue")
# ax.scatter(a[0], a[1], label="a", color = "pink")
# ax.scatter(b[0], b[1], label="b", color = "purple")
ax.set_aspect(aspect='equal')

# ax.legend()

plt.show()