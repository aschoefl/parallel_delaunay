import numpy as np
import ast
from matplotlib import pyplot as plt

import os
import time

N = 10

file_dir = '/home/ams/Studium/Parallel Computations/project/parallel_delaunay'

# submit job
os.chdir(file_dir)
os.system('make')
time.sleep(3)
os.system('./delauney')
time.sleep(2)

myfile = open(file_dir+'/'+'points.txt', 'rt')
data = myfile.read()
pnts = np.array(ast.literal_eval(data))
# pnts = [tuple(int(v) for v in a.strip("()").split(",")) for a in s.split(') (')]
# print(pnts)

plt.scatter(pnts[:,0]+0.5, pnts[:,1]+0.5)
plt.grid()
plt.xticks([ i for i in range(0,N)])
plt.yticks([ i for i in range(0,N)])


plt.show()