from sklearn.exceptions import NonBLASDotWarning
import gaussian
import numpy as np

class Hagedorn(gaussian.Gaussian):

    def __init__(self, q, p, Q, P, eps, c=np.array([1]), s=0, state=None):
        gaussian.Gaussian.__init__(self, q, p, Q, P, eps, s, state)
        
        # set coefficients of Hagedorn wavepackets
        self.c = c

    def get_polynomial(self, x, k):
        y = np.ones(x.shape, dtype=np.complex128)
        z = np.ones(x.shape, dtype=np.complex128)
        temp = np.ones(x.shape, dtype=np.complex128)
        for ix in range(k):
            temp = (np.sqrt(2/self.eps)/self.Q*(x-self.q)*y - np.conj(self.Q)/self.Q*np.sqrt(ix)*z)/np.sqrt(ix + 1)
            z = y
            y = temp

        return y

    def get_hagedorn(self, x, k):

        return self.get_polynomial(x, k) * gaussian.Gaussian.psi(self,x)

    def psi(self, x, n = None):
        # could you use the map function

        # sum over wavepackets
        sum = np.zeros(x.shape, dtype=np.complex128)
        if n == None:
            for k, c in enumerate(self.c):
        #for k, c in [(0, self.c[0])]:
                sum += c*self.get_hagedorn(x, k)
        else:
            for k, c in enumerate(self.c[:n]):
                sum += c*self.get_hagedorn(x, k)
        return np.exp(-1j/self.eps * self.s) * sum

