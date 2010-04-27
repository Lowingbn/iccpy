import numpy.random
import numpy as np

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
    pos, r = sample_density(halo, n, method)
    n = pos.shape[0]
    
    psi = -halo.potential(r)
        
    v_scale = np.sqrt(1 + (r/halo.anisotropy_radius)**2)
    v_perpendicular_max = np.sqrt(2*psi)/v_scale
    prl = pos/r.reshape(len(r), 1)
    vel = np.empty([n,3])
        
    for i in range(n):
        fQ_max = halo.distribution_function_Q(psi[i])
        while True:
            v_i = v_perpendicular_max[i]*numpy.random.rand()**(1.0/3.0)*random_direction()[0]
            v_parallel = np.dot(v_i, prl[i]) 
            v_i = v_i + v_parallel * (v_scale[i] - 1) * prl[i]
                
            L2 = np.square((pos[i,1]*v_i[2] - pos[i,2]*v_i[1], pos[i,2]*v_i[0] - pos[i,0]*v_i[2], pos[i,0]*v_i[1] - pos[i,1]*v_i[0])).sum()
            v2 = np.square(v_i).sum()
            Q = psi[i] - 0.5*v2
            if halo.anisotropy_radius!=np.inf:
                Q -= 0.5*L2/halo.anisotropy_radius**2 
            
            fQ = halo.distribution_function_Q(Q)
    
            if numpy.random.rand()*fQ_max <= fQ: break
        
        vel[i,:] = v_i
            
    return pos, vel

if __name__=="__main__":
    import iccpy.constants
    import iccpy.halo.hernquist
    import iccpy.halo.analysis
    import pylab as pl
    
    #Use Gadget units
    iccpy.constants.set_units('GALACTIC')
    #Create a hernquist profile
    hernquist = iccpy.halo.hernquist.Hernquist(100, 0.03, 0.04)
    
    #Sample the profile
    pos, vel = sample(hernquist, 10000, method='random')
    particle_mass=hernquist.total_mass/float(pos.shape[0])
    
    #Work out density profile of samples
    density, r = iccpy.halo.analysis.density_profile(pos, particle_mass=particle_mass)
    pl.loglog(r, density, label='sampled')
    
    #Compare against actual density profile
    density = hernquist.density(r)
    pl.loglog(r, density, label='analytic')
    
    pl.ylabel(r'$\rho$')
    pl.xlabel(r'$r$ [ $\mathrm{Mpc}$ ]')
    
    pl.legend()
    pl.figure()
    
    #Work out velocity anisotropy profile
    beta, r = iccpy.halo.analysis.anisotropy_profile(pos, vel, n_bins=20)
    pl.semilogx(r, beta, label='sampled')
    
    #Compare against desired profile
    beta = [ x**2 / (x**2 + hernquist.anisotropy_radius**2) for x in r ]
    pl.semilogx(r, beta, label='analytic')
    
    pl.ylabel(r'$\beta$')
    pl.xlabel(r'$r$ [ $\mathrm{Mpc}$ ]')
    pl.legend()
    pl.ylim(-1.5,1.5)
    pl.show()