import numpy as np
import ast
from matplotlib import pyplot as plt
from numpy.lib.financial import _ppmt_dispatcher
from scipy.spatial import Delaunay
import sys


# v = (0.0486388, 0.0646258)
# center = (0.0314435, 0.0182392)
# points = [ (0, 1), (0, 0), (1, 0), ]
# voronoi = [ (-0.26902, 0.5), (0.5, -0.825753), (0.525445, 0.525445), ]
# to remove : 0,2

v =  (0.30402, 0.427112)
center = (0.269364, 0.438579)
a = (0.286692, 0.432846)
b = (0.280959, 0.415517)
points = [ (0.311819, 0.457669), (0.294037, 0.469799), (0.226632, 0.521761), (0.20404, 0.45074), (0.245476, 0.365896), (0.247831, 0.362722), (0.314423, 0.393222), ]
voronoi =  [ (0.291262, 0.446633), (0.248499, 0.480427), (0.242817, 0.477508), (0.230454, 0.4111), (0.2854, 0.393042), (0.272555, 0.396688), (0.301009, 0.424956), ]


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

if( (a[0]-b[0]) != 0):
    h = lambda x: (a[1]-b[1])/(a[0]-b[0])*(x-b[0])+b[1]
    ax.plot(np.linspace(0,0.5), h(np.linspace(0,0.5)), label = "$h_v$")
# f = lambda x: (center[1]-v[1])/(center[0]-v[0])*(x-v[0])+v[1]



# ax.plot(np.linspace(0,0.5), f(np.linspace(0,0.5)))


ax.plot(voronoi[:,0], voronoi[:,1], '-o', label="Q", color="orange")
ax.scatter(center[0], center[1], label="c", color="red")
ax.scatter(v[0], v[1], label="v", color="blue")
# ax.scatter(a[0], a[1], label="a", color = "pink")
# ax.scatter(b[0], b[1], label="b", color = "purple")
ax.set_aspect(aspect='equal')
plt.plot(points[:,0], points[:,1], '-o', label="P", color="lime")
ax.legend()

plt.show()