import numpy as np 
import sys

sys. stdout = open("test_sample_py.txt", "w")
K = 200
eps = 0.1
mean = 0
std = np.sqrt(eps/2)

print("npart, rmse")
for i in range(1,8):
    npart = 10**i
    rmse = 0
    for k in range(K):
        sum_ = np.sum(np.random.normal(0,std, npart))/npart
        rmse += sum_**2
    print(npart, np.sqrt(rmse/K))
sys. stdout. close()



