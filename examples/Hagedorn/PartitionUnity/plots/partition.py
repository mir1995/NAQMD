import numpy as np 

def partition_f(x, mean, delta):

    diff = np.abs(x - mean)

    mask = diff < delta / 2

    diff[~mask] = 0
    a = (delta / 2) / (- diff[mask] + delta / 2)
    b = (delta / 2) / diff[mask]
    diff[mask] = np.exp(-a) / (np.exp(-a) + np.exp(-b))

    return diff
