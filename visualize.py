import numpy as np
import ast
from matplotlib import pyplot as plt

import os
import time

N = 20

file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay'

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

myfile = open(file_dir+'/'+'polyPoints0.txt', 'rt')
data = myfile.read()
pnts = np.array(ast.literal_eval(data))
plt.scatter(pnts[:,0], pnts[:,1])

myfile = open(file_dir+'/'+'veronoiPoints0.txt', 'rt')
data = myfile.read()
pnts = np.array(ast.literal_eval(data))

plt.scatter(pnts[:,0], pnts[:,1])
plt.grid()
plt.xticks([ i/float(N) for i in range(0,N)])
plt.yticks([ i/float(N) for i in range(0,N)])


plt.show()