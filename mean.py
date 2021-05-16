from matplotlib import pyplot as plt
import numpy as np
from numpy.core.records import array
from numpy.testing._private.utils import print_assert_equal 

no_pnts = [368277]
max_time = [65125837/1000]
R = [1]

no_pnts += [92088+92042+86677+91784]
max_time += [24023419/1000]

R += [4]

no_pnts += [40888+40912+40937+40968+40790+41009+40911+40977+40843]
max_time += [10851959/1000]
R += [9]


no_pnts += [23075+22998+23097+23066+22933+23001+23012+\
    23042+23011+22947+23001 +23013+23020+22952+23039+22923]
max_time += [5709119/1000]
R += [16]

no_pnts = np.array(no_pnts)
max_time = np.array(max_time)
R = np.array(R)

print("R = ", R)
print("no_pnts = ", no_pnts)
print("max_time = ", max_time)
print("time/point = ", max_time/no_pnts)
print("speedup = ", (max_time/no_pnts)[0]/(max_time/no_pnts))

plt.plot(R, max_time/no_pnts, '-o')
plt.xlabel("amount of processors")
plt.ylabel("computation time per point [ms]")

plt.xscale('log')

plt.show()



