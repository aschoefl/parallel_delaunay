import numpy as np
import ast
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay

import os
import time
import sys

N = 10
P = 1


file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay/outdata/'

## submit job
# os.chdir(file_dir)
# os.system('make')
# time.sleep(3)
# os.system('mpiexec -n 4 delauney')
# time.sleep(2)

# for k in range(0,P*P):
#     myfile = open(file_dir+'points'+str(k)+'.txt', 'rt')
#     data = myfile.read()
#     if k == 0:
#         points = np.array(ast.literal_eval(data))
#     else:
#         points = np.concatenate([points, np.array(ast.literal_eval(data))])

#     myfile.close()

scale = 1
fig = plt.figure()
ax = fig.add_subplot(111)
# plt.figure()
# ax.scatter(points[:,0], points[:,1], label="pnts", color="grey")
plt.grid()
ax.set_aspect(aspect='equal')
plt.xticks([ scale*i/float(N) for i in range(0,N)])
plt.yticks([ scale*i/float(N) for i in range(0,N)])
# plt.legend()
# plt.show()

ind_points = 0
for i in range(0,N):
    for j in range(0,N):
        ind_delauney = i*1000+j*10

        # plt.figure(ind_delauney)
        # plt.scatter(points[:,0], points[:,1], label="pnts", color="grey")

        myfile = open(file_dir+'polyPoints'+str(ind_delauney)+'.txt', 'rt')
        data = myfile.read()
        pnts = np.array(ast.literal_eval(data))*scale
        center = pnts[0,:]
        for k in range(len(pnts[1:,0])):
            ax.plot([pnts[0,0], pnts[k,0]], [pnts[0,1], pnts[k,1]],'-o', label="poly", color="black")
            # plt.scatter(pnts[0,0], pnts[0,1], label="center", color="red")


        # myfile = open(file_dir+'voronoiPoints'+str(ind_delauney)+'.txt', 'rt')
        # data = myfile.read()
        # pnts = np.array(ast.literal_eval(data))
        # plt.plot(pnts[:,0], pnts[:,1], label="voronoi", color = "orange")


        # plt.grid()
        # plt.xticks([ i/float(N) for i in range(0,N)])
        # plt.yticks([ i/float(N) for i in range(0,N)])
        # plt.legend()
        # break
    # plt.show()

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