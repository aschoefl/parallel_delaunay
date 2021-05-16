import numpy as np
import ast
from matplotlib import pyplot as plt
from scipy.spatial import Delaunay

N = 8
P = 2

file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay/outdata/'


fig = plt.figure()
ax = fig.add_subplot(111)
points = []
ind_points = 0
for i in range(0,N):
    for j in range(0,N):
        ind_delauney = i*1000+j*10

        myfile = open(file_dir+'polyPoints'+str(ind_delauney)+'.txt', 'rt')
        data = myfile.read()
        pnts = np.array(ast.literal_eval(data))
        center = pnts[0,:]
        for k in range(len(pnts[1:,0])):
            ''' delete ghost points '''
            if (pnts[k,0] in [0,1] and pnts[k,1] in [0,1]):
                continue
            ax.plot([pnts[0,0], pnts[k,0]], [pnts[0,1], pnts[k,1]],'-o', label="poly", color="black")
            # plt.scatter(pnts[0,0], pnts[0,1], label="center", color="red")
        
        points += [tuple(center)]


plt.grid()
ax.set_aspect(aspect='equal')
plt.xticks([ i/float(N) for i in range(0,N)])
plt.yticks([ i/float(N) for i in range(0,N)])
ax.set_xlim([0,1])
ax.set_ylim([0,1])
plt.show()

points = np.array(points)
if True:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tri = Delaunay(points)
    ax.triplot(points[:,0], points[:,1], tri.simplices.copy(), color = 'black')
    ax.scatter(points[:,0], points[:,1], color='black')
    plt.grid()
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_aspect(aspect='equal')
    plt.xticks([ i/float(N) for i in range(0,N)])
    plt.yticks([ i/float(N) for i in range(0,N)])
    plt.show()