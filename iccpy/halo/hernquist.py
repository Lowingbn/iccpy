import numpy as np
import iccpy.constants as constants

class Hernquist:
    def __init__(self, total_mass, scale_length, r_a=np.inf):
        self.total_mass = total_mass
        self.scale_length = scale_length
        self.anisotropy_radius = r_a
        
    def density(self, r):
        return self.total_mass*self.scale_length/(2*np.pi*r*(r+self.scale_length)**3)
    
    def cumlative_mass(self, r):
        return self.total_mass*r**2 /(r+self.scale_length)**2
    
    def radius(self, mass):
        m = np.sqrt(mass / self.total_mass)
        return self.scale_length * m / (1 - m)
    
    def potential(self, r):
        return -constants.G*self.total_mass/(r+self.scale_length)
    
    def distribution_function_Q(self, Q):
        M = self.total_mass
        a = self.scale_length
        
        vg = np.sqrt(constants.G*M/a)
        q = np.sqrt(a*Q/(constants.G*M))
    
        ans  = M * (3*np.arcsin(q) + q*(1-q**2)**0.5 * (1-2*q**2) * (8*q**4-8*q**2-3))
        ans /= (8*np.sqrt(2)*np.pi**3*a**3*vg**3) * (1-q**2)**2.5
        
        if self.anisotropy_radius!=np.inf:
            ans2 = M * a**2 * q * (1-2*q**2) 
            ans2 /= np.sqrt(2) * np.pi**3*a**3*vg**3 * self.anisotropy_radius**2
            ans += ans2
    
        return ans
    
    def distribution_function(self, pos, vel):
        r = np.square(pos).sum()
        psi = -self.potential(r)
        L2 = np.square((pos[1]*vel[2] - pos[2]*vel[1], pos[2]*vel[0] - pos[0]*vel[2], pos[0]*vel[1] - pos[1]*vel[0])).sum()
        v2 = np.square(vel).sum()
        Q = psi - 0.5*v2 - 0.5*L2/self.anisotropy_radius**2 
        
        return self.distribution_function(Q)