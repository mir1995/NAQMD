import numpy as np

class Gaussian():

    def __init__(self, q, p, Q, P, eps, s=0, state=None):

        self.q = q
        self.p = p
        self.Q = Q
        self.P = P
        self.eps = eps
        self.q = q
        self.s = s 
        self.state = state

    def psi(self, x):

        return (np.pi * self.eps)**(-0.25) * self.Q**(-0.5) * np.exp(1j / 2 / self.eps * self.P / self.Q * (x - self.q)**2 + 1j*self.p*(x-self.q))

