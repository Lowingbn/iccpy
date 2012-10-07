import numpy as np
import iccpy.constants as constants

def density(r, rho_s, r_s):
    x = r/r_s
    return rho_s/(x*(1+x)**2)

class NFW:
    def __init__(self, rho_s, scale_length):
        self.rho_s = rho_s
        self.scale_length = scale_length
        
    def density(self, r):
        x = r/self.scale_length
        return self.rho_s/(x*(1+x)**2)
