import numpy as np
import pylab as pl
import scipy.optimize

def density_profile(pos, particle_mass=1, n_bins=20, centre=None, r_min=None, r_max=None):
    """Calculates the spherically averaged density profile of a set of particles
    
    Keyword arguments:
    pos           -- n x 3 array of positions
    particle_mass -- scalar or array (default 1.0)
    n_bins        -- number of bins to use (default 20)
    centre        -- if set then used as the origin (default None)
    r_min         -- the start radius for the profile (default None)
    r_max         -- the outer radius for the profile (default None)    
    """
    if centre is None:
        r_pts = np.sqrt(np.square(pos).sum(1))
    else:
        r_pts = np.sqrt(np.square(pos-centre).sum(1))
    
    #Find the minimum r that is non-zero   
    if r_min is None:
        r_min = np.sort(r_pts)[1]
    if r_max is None:
        r_max = r_pts.max()
        
    r_bins = r_min * np.power(r_max/r_min, np.arange(n_bins)/(n_bins-1.0))
    
    if isinstance(particle_mass, np.ndarray):
        #Each particle has it's own mass
        bin_mass, r_bins = np.histogram(r_pts, r_bins, weights=particle_mass)
        bin_vol = r_to_dV(r_bins)[1:]
        bin_dens = bin_mass / bin_vol
    else:
        #All particles have the same mass
        bin_mass, r_bins = np.histogram(r_pts, r_bins)
        bin_vol = r_to_dV(r_bins)[1:]
        bin_dens = (bin_mass * particle_mass) / bin_vol
    
    return bin_dens, r_bins[1:] - np.diff(r_bins)/2
    
def r_to_dV(rad):
    """ find differential vol from radii """
    vol = 4.0/3.0 * np.pi * rad * rad * rad
 
    dV = np.empty_like(vol)
    dV[0], dV[1:] = vol[0], np.diff(vol)

    return dV 
    
def anisotropy_profile(pos, vel, n_bins=20, centre=None):
    if centre is None:
        rpos = pos    
    else:
        rpos = pos - centre
    
    r_pts = np.sqrt(np.square(rpos).sum(1))
    r_min = np.sort(r_pts)[1]
    r_bins = r_min * np.power(r_pts.max()/r_min, np.arange(n_bins)/(n_bins-1.0))
    
    dir = rpos/r_pts.reshape(r_pts.shape[0], 1)
    
    v_r = (dir*vel).sum(1)
    v2 = np.square(vel).sum(1)
    v_p = np.sqrt(v2 - v_r**2)
    
    bin_nums = np.histogram(r_pts, r_bins)[0]
    #bin_v_r = np.histogram(r_pts, r_bins, weights=v_r)[0]
    bin_v_r2 = np.histogram(r_pts, r_bins, weights=v_r**2)[0]
    #bin_v_p = np.histogram(r_pts, r_bins, weights=v_p)[0]
    bin_v_p2 = np.histogram(r_pts, r_bins, weights=v_p**2)[0]
         
    #mean_v_r = bin_v_r/bin_nums
    sigma_v_r = bin_v_r2/bin_nums #- mean_v_r**2
    
    #mean_v_p = bin_v_p/bin_nums
    sigma_v_p = bin_v_p2/bin_nums #- mean_v_p**2
    
    beta = 1 - sigma_v_p/(2*sigma_v_r)
    
    return beta, r_bins[1:] - np.diff(r_bins)/2
    
def _nfw(params, radii, density):
    x = radii/params[1]
    dfdp0 = 1/(x*(1+x)**2)
    f = params[0]*dfdp0
    
    dxdp1 = -radii/params[1]**2
    dfdx = -params[0]*(3*x+1)/(x**2*(1+x)**3)
    dfdp1 = dxdp1*dfdx
    
    J = np.sum((density - f)**2)
    dJdp0 = 2*np.sum(dfdp0*(f - density))
    dJdp1 = 2*np.sum(dfdp1*(f - density))

    return J, np.array([dJdp0, dJdp1])
    
def fit_nfw(radii, density):
    r_s0 = np.sqrt(density[-1]*radii[-1]**3 / (density[0]*radii[0]))
    rho_s0 = radii[0]*density[0]/r_s0
            
    guess = (rho_s0, r_s0)        
    bounds = ((0, None), (0, None))
    ret = scipy.optimize.fmin_tnc(_nfw, guess, approx_grad=False, bounds=bounds, disp=0, args=(radii, density))
    
    if ret[2]==1 or ret[2]==2:
        return ret[0]
    else:        
        return None
        
def plot_nfw(rho_0, r_s, r_min=None, r_max=None):
    if r_min==None: r_min = 0.01*r_s
    if r_max==None: r_max = 100 * r_s
    
    r = np.logspace(np.log10(r_min), np.log10(r_max), 1000)
    x = r/r_s
    rho = rho_0 / (x*(1+x)**2)
    
    pl.loglog(r, rho)
    
if __name__=="__main__":
    import pylab as pl
    
    r = np.logspace(-2,3,1000)
    rho = 1000/(r*(1+r)**2)
    
    params = fit_nfw(r, rho)
    plot_nfw(*params)
    pl.show()
    
