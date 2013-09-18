import numpy as np

def calculate_halo_velocity(pos, vel, centre, radius, masses=None):
    r = np.sqrt(np.sum((pos-centre)**2, axis=1))
    idx = np.where(r<radius)[0]

    if masses is None:
        return np.average(vel[idx], axis=0)
    else:
        return np.average(vel[idx], axis=0, weights=masses[idx])


