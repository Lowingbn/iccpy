import numpy as np

def density_profile(pos, particle_mass=1, n_bins=20, centre=None):
    if centre is None:
        r_pts = np.sqrt(np.square(pos).sum(1))
    else:
        r_pts = np.sqrt(np.square(pos-centre).sum(1))
   
    r_min = np.sort(r_pts)[1]
    r_bins = r_min * np.power(r_pts.max()/r_min, np.arange(n_bins)/(n_bins-1.0))

    bin_mass, r_bins = np.histogram(r_pts, r_bins)
    bin_vol = r_to_dV(r_bins)[1:]
    bin_dens = (bin_mass * particle_mass) / bin_vol
    
    return bin_dens, r_bins[1:] - np.diff(r_bins)/2, r_bins
    
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

def r_to_dV(rad):
    """ find differential vol from radii """
    vol = 4.0/3.0 * np.pi * rad * rad * rad
 
    dV = np.empty_like(vol)
    dV[0], dV[1:] = vol[0], np.diff(vol)

    return dV
