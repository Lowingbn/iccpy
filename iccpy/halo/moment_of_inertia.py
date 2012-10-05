import numpy as np
import scipy.linalg

def moi(pos, particle_mass, max_radius, reduced=False, centre=None):
    """Calculates the moment of inertia of a set of particles
    
    Keyword arguments:
    pos           -- n x 3 array of positions
    particle_mass -- either a scalar or a 1D length n array of individual particle Masses
    max_radius    -- maximum radius of particles for which in include in the MoI
    reduced       -- whether to use the reduced or the full MoI (default False)
    centre        -- if set then used as the origin (default None)
    """

    if centre is None: pos_i = pos
    else: pos_i = pos - centre

    r_pts = np.sqrt(np.square(pos_i).sum(1))
    idx = np.where(r_pts<=max_radius)

    #First
    moi = np.zeros([3,3])
    
    if reduced:
        if isinstance(particle_mass, np.ndarray):
            for i in idx:
                moi += particle_mass[i] * np.outer(pos_i[i], pos_i[i]) / np.dot(pos_i[i], pos_i[i])
        else:
            for i in idx:
                moi += np.outer(pos_i[i], pos_i[i]) / np.dot(pos_i[i], pos_i[i])
            moi *= particle_mass
    else:
        if isinstance(particle_mass, np.ndarray):
            for i in idx:
                moi += particle_mass[i] * np.outer(pos_i[i], pos_i[i])
        else:
            for i in idx:
                moi += np.outer(pos_i[i], pos_i[i])
            moi *= particle_mass

    values, orient = scipy.linalg.eigh(moi)

    return orient[:,0]
