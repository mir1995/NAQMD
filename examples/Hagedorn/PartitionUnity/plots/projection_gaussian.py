import numpy as np 


def project(Wave, q):
    """
    Wave: Gaussian incoming wavepacket

    Consider two Gaussians 

    Their product is given by https://en.wikipedia.org/wiki/Gaussian_integral

    """
    A = np.zeros(shape=(q.shape, q.shape))
    
    b = np.zeros(q.shape)

    for i, x in enumerate(q):
        b[i] = 

