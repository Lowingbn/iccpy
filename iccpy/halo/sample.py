import numpy.random
import numpy as np
import scipy.integrate

import iccpy.constants as constants

def _mass_integrand(r, halo):
    if (r==0): return 0
    else: return r*r*halo.density(r)
    
def _potential_integrand(r, halo):
    return r*halo.density(r)

class HaloRealisation:
    def __init__(self, alpha, beta, gamma, rho_s, r_s, r_a, r_v=0):
        print "Initialising halo model..."
        
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        
        self.rho_s = rho_s
        self.r_s = r_s
        self.anisotropy_radius = r_a
        
        self.d = (self.beta - self.gamma) / self.alpha  #tested
        
        if beta<=3:
            self.r_vir = r_v
            self.r_decay = 10*self.r_vir
            self.c = self.r_vir / self.r_s
            self.eplison = -(self.gamma + self.beta * self.c**self.alpha)/(1 + self.c**self.alpha) \
                            + (self.r_vir/self.r_decay)
            self.density_const = self.rho_s / ((self.c**self.gamma)*((1+self.c**self.alpha)**self.d))
            r_max = self.r_vir+1000*self.r_decay
        else:
            r_max = 1e5*self.r_s
            
        r_min = self.r_s * 1e-6
        self.r_table = np.logspace(np.log10(r_min),np.log10(r_max),1e4)
        
        print "Tabulating mass..."
        self.mass_table = np.array([ self.cumlative_mass_actual(r) for r in self.r_table ])
        
        print "Tabulating potential..."
        self.potential_table = np.array([ self.potential_actual(r) for r in self.r_table ])
        
        print "Tabulating f(Q)... warning this may take some time"
        self.logQ_table = np.linspace(np.log10(-self.potential_table[-1000]), np.log10(-self.potential_table[0]), 10000, endpoint=True)        
        self.logf_table = np.log10( [ self.dist_func_actual(10**logQ) for logQ in self.logQ_table ])
        
        #Choose an values where we know that f(Q) is still correct
        self.r_max = self.r_table[-1000]
        self.total_mass = self.mass_table[-1000]
        
        print "Initialisation Complete"
                        
    def density(self, r): #tested
        #if r>0.04: return 0
        #print r
        if self.beta<=3 and r>self.r_vir:
            return self.density_const * ((r/self.r_vir)**self.eplison) * np.exp(-(r-self.r_vir)/self.r_decay)
        else:
            x = r/self.r_s
            return self.rho_s / ((x**self.gamma) * ((1 + x**self.alpha)**self.d))
    
    def cumlative_mass_actual(self, r):  #tested
        return 4*np.pi*scipy.integrate.quad(_mass_integrand, 0, r, (self))[0]
    
    def cumlative_mass(self, r):
        return np.interp(r, self.r_table, self.mass_table)
    
    def radius(self, mass):
        return np.interp(mass, self.mass_table, self.r_table)

    def potential_actual(self, r):  #tested
        ans = constants.G * self.cumlative_mass_actual(r) /r
        points = np.logspace(np.log10(r), 10+np.log10(self.r_s), 10, False)[1:]
        ans += 4*np.pi*constants.G*scipy.integrate.quad(_potential_integrand, r, 1e10*self.r_s, (self), points=points)[0]
        return -ans
    
    def potential(self, r):  #tested
        return np.interp(r, self.r_table, self.potential_table)

    def dpsidr(self, r):  #tested
        return -constants.G * self.cumlative_mass(r) / r**2
    
    def dpsidr_actual(self, r):  #tested
        return -constants.G * self.cumlative_mass_actual(r) / r**2
    
    def d2psidr2(self, r): #tested
        return -4*np.pi*constants.G*self.density(r) - 2*self.dpsidr(r)/r
    
    def rhoQ(self, r): #tested
        return self.density(r) * (1 + (r/self.anisotropy_radius)**2)
    
    def drhoQdr(self, r): #tested
        return self.drhodr(r)*(1+(r/self.anisotropy_radius)**2) + 2*r*self.density(r)/self.anisotropy_radius**2
    
    def d2rhoQdr2(self, r): #tested
        return self.d2rhodr2(r)*(1+(r/self.anisotropy_radius)**2) + 4*r*self.drhodr(r)/self.anisotropy_radius**2 + 2*self.density(r)/self.anisotropy_radius**2
    
    def drhoQdpsi(self, r): #tested
        return self.drhoQdr(r)/self.dpsidr(r)
    
    def d2rhoQdpsi2(self, r): #tested        
        dsdr = self.dpsidr(r)
        return (self.d2rhoQdr2(r) / dsdr**2) - (self.d2psidr2(r) * self.drhoQdr(r) / dsdr**3)
            
    def drhodr(self, r): #tested
        if self.beta<=3 and r>self.r_vir:
            e = self.eplison
            expr = np.exp(-(r-self.r_vir)/self.r_decay)
            ans = e/r - 1/self.r_decay
            return self.density_const * (r/self.r_vir)**e * ans * expr
        else:
            a = self.rho_s
            b = self.gamma
            c = self.alpha
            d = self.d
            x = r/self.r_s
            
            ans =  -a*b*(x**-(b+1))*(1+x**c)**-d -a*c*d*(x**-b)*(1+x**c)**-(d+1)
            return ans/self.r_s
        
    def d2rhodr2(self, r): #tested
        if self.beta<=3 and r>self.r_vir:
            e = self.eplison
            er = np.exp(-(r-self.r_vir)/self.r_decay)
            
            ans =  e*(e-1)*r**(e-2)*er / self.r_vir**e
            ans -= 2*e*r**(e-1)*er/(self.r_vir**e * self.r_decay)
            ans += (r/self.r_vir)**e * er / self.r_decay**2

            return ans*self.density_const
        else:
            a = self.rho_s
            b = self.gamma
            c = self.alpha
            d = self.d
            x = r/self.r_s
            xc = (1+x**c)
            
            ans =  a*b*(b+1)*(x**-(b+2))*xc**-d
            ans += a*b*c*d*(x**-(b-c+2))*xc**-(d+1)
            ans += a*c*d*(b-c+1)*x**-(b-c+2)*xc**-(d+1)
            ans += a*c*c*d*(d+1)*x**-(b-2*c+2)*xc**-(d+2)
            
            return ans/self.r_s**2
    
    def r_psi(self, psi):
        return np.interp(-psi, self.potential_table, self.r_table)
    
    def dist_func_integrand(self, x, Q):
        r = self.r_psi(Q *np.sin(x)**2)
        return self.d2rhoQdpsi2(r)# * np.sin(x)
    
    def dist_func_actual(self, Q):
        x_min = np.arcsin(np.sqrt(-self.potential_table[-1]/Q))
        fQ = 2*np.sqrt(Q)*scipy.integrate.quad(self.dist_func_integrand, x_min, np.pi/2, (Q), epsabs=0, epsrel=1e-5, weight='sin', wvar=1)[0]/(np.sqrt(8)*np.pi**2)
        print Q, fQ
        return fQ
    
    def distribution_function_Q_max(self, Q):
        logQ = np.log10(Q)
        idx = np.searchsorted(self.logQ_table, logQ)
        return 10**np.array([ np.max(self.logf_table[:i]) for i in idx])
    
    def distribution_function_Q(self, Q):
        logQ = np.log10(Q)
        logf = np.interp(logQ, self.logQ_table, self.logf_table, left=0)
        return 10**logf
      
def random_direction(n=1):
    phi = numpy.random.rand(n) * 2.0 * np.pi

    costheta = numpy.random.rand(n) * 2 - 1;  
    sintheta = np.sqrt(1 - costheta * costheta)
    return np.column_stack((np.sin(phi)*sintheta, np.cos(phi)*sintheta, costheta))

def sample_density(halo, n, method='random'):
    if method=='random':
        mass = halo.total_mass * numpy.random.rand(n)
        r = halo.radius(mass)
        pos = random_direction(n)*r.reshape(len(r), 1)
        
        return pos, r
    elif method=='grid':
        n_x = np.ceil((6*n/np.pi)**(1.0/3.0))
        if n_x%2==1: n_x+=1
        pos = 2*((np.mgrid[:n_x, :n_x, :n_x] + 0.5)/ n_x) - 1
        pos = pos.reshape((3, n_x**3)).swapaxes(0,1)
            
        r2 = np.square(pos).sum(1)
        idx = np.flatnonzero(r2<=1.0)

        pos = pos[idx,:]
        r_old = np.sqrt(r2[idx])
            
        mass = halo.total_mass * r_old**3
        r = halo.radius(mass)
        pos = pos*(r/r_old).reshape(len(r),1)
        n = pos.shape[0]
        
        return pos, r
    else:
        return None
    
def sample(halo, n, method='random'):
    print "Sampling halo"
    pos, r = sample_density(halo, n, method)
    n = pos.shape[0]
    
    print "Density sampling complete, %d points picked" % n
    
    psi = -halo.potential(r)
        
    v_scale = np.sqrt(1 + (r/halo.anisotropy_radius)**2)
    v_perpendicular_max = np.sqrt(2*psi)/v_scale
    prl = pos/r.reshape(len(r), 1)
    prl = prl.swapaxes(0,1)
    vel = np.empty([n,3])

    fQ_max = halo.distribution_function_Q_max(psi)        

    idx = np.arange(n, dtype=np.int32)
    rejections = 0

    while idx.size > 0:
        v_i = v_perpendicular_max[idx]*numpy.random.random(idx.size)**(1.0/3.0)*random_direction(idx.size).swapaxes(0,1)
        
        v_parallel = (prl[:,idx] * v_i).sum(0) 
        
        v_i = v_i + v_parallel * (v_scale[idx] - 1) * prl[:,idx]
                
        L2 = np.square((pos[idx,1]*v_i[2] - pos[idx,2]*v_i[1], pos[idx,2]*v_i[0] - pos[idx,0]*v_i[2], pos[idx,0]*v_i[1] - pos[idx,1]*v_i[0])).sum(0)

        v2 = np.square(v_i).sum(0)
        Q = psi[idx] - 0.5*v2

        if hasattr(halo, 'anisotropy_radius') and halo.anisotropy_radius!=np.inf:
            Q -= 0.5*L2/halo.anisotropy_radius**2 
        
        fQ = halo.distribution_function_Q(Q)
        
        reject = numpy.random.random(idx.size) * fQ_max[idx] > fQ
        keep = np.flatnonzero(np.logical_not(reject))

        vel[idx[keep]] = v_i[:,keep].swapaxes(0,1)        
        idx = idx[np.flatnonzero(reject)]

        rejections += idx.size
        
            
    print 'Velocity sampling done, total rejected', rejections

    return pos, vel

def sample_hernquist(num_particles, mass, scale_length, anisotropy_radius=np.inf):
    import iccpy.halo.hernquist
    
    #Create a hernquist profile
    hernquist = iccpy.halo.hernquist.Hernquist(mass, scale_length, anisotropy_radius)
    return sample(hernquist, num_particles)
