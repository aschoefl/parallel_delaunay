import numpy as np
import ast
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay

import os
import time

N = 4

file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay/outdata'

## submit job
# os.chdir(file_dir)
# os.system('make')
# time.sleep(3)
# os.system('mpiexec -n 4 delauney')
# time.sleep(2)

# for k in range(0,4):
#     myfile = open(file_dir+'/'+'points'+str(k)+'.txt', 'rt')
#     data = myfile.read()
#     if k == 0:
#         pnts = np.array(ast.literal_eval(data))
#     else:
#         pnts = np.concatenate([pnts, np.array(ast.literal_eval(data))])

#     myfile.close()


for k in range(0,4):
    add = N
    plt.figure(k)
    # myfile = open(file_dir+'/'+'points'+str(k)+'.txt', 'rt')
    # data = myfile.read()
    # pnts = np.array(ast.literal_eval(data))
    # plt.scatter(pnts[:,0], pnts[:,1], label="pnts", color="grey")
    # points = pnts



    try:
        myfile = open(file_dir+'/'+'polyPoints'+str(k+add)+'.txt', 'rt')
    except:
        continue
    data = myfile.read()
    pnts = np.array(ast.literal_eval(data))
    center = pnts[0,:]
    plt.scatter(pnts[1:,0], pnts[1:,1], label="poly", color="lime")
    plt.scatter(pnts[0,0], pnts[0,1], label="center", color="red")
    # if k==3: plt.scatter(0.375, 0.875)
    # if k==1: 
    #     plt.scatter(0.875, 0.875)
        # plt.scatter((0.875+pnts[0,0])/2, (0.875+pnts[0,1])/2)
    # v = np.array([(0.5, 0.75), (0, -0.125), (0.375, 0.625), (0.625, 0.375), (0.875, 0.625)])
    # plt.plot(v[:,0], v[:,1])

    myfile = open(file_dir+'/'+'voronoiPoints'+str(k+add)+'.txt', 'rt')
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