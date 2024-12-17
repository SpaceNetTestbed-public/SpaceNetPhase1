import numpy as np

class DoTD_History:

    def __init__(self, M, T):
        self.M = M
        self.T = T
        self.PI_i_t = np.zeros(((M, self.T+1)))
        self.phi_i_j_t = np.zeros(((M, self.M, self.T+1)))
        self.S_max = np.zeros((self.T+1))
        self.L_max = np.zeros((self.T+1))
        self.t = 0
        self.w1 = 0
        self.w2  = 0
        self.w3 = 1-(self.w1 + self.w2)

    def step(self):
        self.t += 1
