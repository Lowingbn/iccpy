import numpy as np
import scipy.spatial
import iccpy.gadget.kernels

class DensityField(object):
    def __init__(self, pos, mass):
        self.pos = pos
        self.mass = mass
        self.tree = scipy.spatial.cKDTree(pos)

    def density(self, pos, num_neighbours=64):
        if len(pos.shape)==1:
            p_pos = pos.reshape(1, pos.shape[0])
        else:
            p_pos = pos

        r, idx = self.tree.query(p_pos, k=num_neighbours, eps=0.0)

        weights = iccpy.gadget.kernels.sph_kernel(r, r[:,-1][:,np.newaxis])
        density = np.sum(weights * self.mass[idx], axis=1)

        return density

if __name__=="__main__":
    import iccpy.halo.sample
    import iccpy.halo.hernquist
    import iccpy.constants
    iccpy.constants.set_units("GALACTIC")

    num_particles = 1000000
    halo = iccpy.halo.hernquist.Hernquist(10, 0.01)
    p_mass = np.ones(num_particles) * halo.total_mass/num_particles
    p_pos = iccpy.halo.sample.sample_density(halo, num_particles)[0]
    
    df = DensityField(p_pos, p_mass)
    
    #test_pos = np.array([0.01, 0, 0])
    #res = []
    #r = np.logspace(-3, 
    
    test_pos = p_pos[::10000,:]
    den = df.density(test_pos, 64)
    r = np.sqrt(np.sum(test_pos**2, axis=1))
    
    import matplotlib.pyplot as pl
    pl.loglog(r, den, 'x')
    r = np.logspace(-3, 0, 100)
    pl.plot(r, halo.density(r))
    pl.show()

    #for i in range(1000):
    #    test_pos = np.array([0.01, 0, 0])
    #    den = df.density(test_pos, 32)
    #    print i, den
    #    res.append(den)

    #print halo.density(0.01)

    #res = np.array(res)
    #import matplotlib.pyplot as pl
    #pl.plot(range(5, 1000), res)
    #pl.show()


