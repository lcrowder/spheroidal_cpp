import numpy as np
import matplotlib.pyplot as plt
from scipy import io

# Load the data
# DL_coincident_16_std=np.loadtxt('../build/examples/data/std/DL_coincident_16_std.csv',dtype=np.complex128)

DL_coincident_16_std = np.genfromtxt('../build/examples/data/std/DL_coincident_16_std.csv', 
                                     converters={0: lambda x: x.replace('i','j')},
                                     dtype=str)

DL_coincident_16_std=io.loadmat('../build/examples/data/std/DL_coincident_16_std.csv')
DL_far_16_std=io.loadmat('../build/examples/data/std/DL_far_16_std.csv')
PointChargePotential = io.loadmat('../build/examples/data/PointChargePotential.csv')

DL_coincident = []
DL_far=[]
Solution_near=[]
for p in range(2,17):
    DL_coincident.append(io.loadmat('../build/examples/data/DL_coincident_'+str(p)+'.csv'))
    DL_far.append(io.loadmat('../build/examples/data/DL_far_'+str(p)+'.csv'))
    Solution_near.append(io.loadmat('../build/examples/data/Solution_near_'+str(p)+'.csv'))


