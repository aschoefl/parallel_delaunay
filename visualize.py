import numpy as np
import ast
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay

import os
import time
import sys

N = 1

file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay/outdata/'

## submit job
# os.chdir(file_dir)
# os.system('make')
# time.sleep(3)
# os.system('mpiexec -n 4 delauney')
# time.sleep(2)

for k in range(0,4):
    myfile = open(file_dir+'points'+str(k)+'.txt', 'rt')
    data = myfile.read()
    if k == 0:
        points = np.array(ast.literal_eval(data))
    else:
        points = np.concatenate([points, np.array(ast.literal_eval(data))])

    myfile.close()

plt.figure()
plt.scatter(points[:,0], points[:,1], label="pnts", color="grey")
plt.grid()
plt.xticks([ i/float(N) for i in range(0,N)])
plt.yticks([ i/float(N) for i in range(0,N)])
plt.legend()
plt.show()

ind_points = 3
ind_delauney = 2000
# for k in range(0,4):
if False:
    plt.figure(ind_delauney)

    myfile = open(file_dir+'/'+'polyPoints'+str(ind_delauney)+'.txt', 'rt')
    data = myfile.read()
    pnts = np.array(ast.literal_eval(data))
    center = pnts[0,:]
    plt.scatter(pnts[1:,0], pnts[1:,1], label="poly", color="lime")
    plt.scatter(pnts[0,0], pnts[0,1], label="center", color="red")


    myfile = open(file_dir+'/'+'voronoiPoints'+str(ind_delauney)+'.txt', 'rt')
    data = myfile.read()
    pnts = np.array(ast.literal_eval(data))
    plt.plot(pnts[:,0], pnts[:,1], label="voronoi", color = "orange")


    plt.grid()
    plt.xticks([ i/float(N) for i in range(0,N)])
    plt.yticks([ i/float(N) for i in range(0,N)])
    plt.legend()
    plt.show()

if False:
    plt.figure()
    tri = Delaunay(points)
    plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    plt.scatter(points[:,0], points[:,1])
    plt.scatter(center[0], center[1], label="center", color="red")
    print(center)
    plt.grid()
    plt.xticks([ i/float(N) for i in range(0,N)])
    plt.yticks([ i/float(N) for i in range(0,N)])

    plt.legend()
    plt.show()