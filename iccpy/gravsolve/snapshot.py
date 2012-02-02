import numpy as np
import h5py

def write_gravsolve_file(filename, num_particles, pos, vel, mass, id=None, time=None):
    f = h5py.File('%s.hdf5' % filename, 'w')
    
    num_particles = np.array(num_particles)
    
    if id is None:
        ids = np.arange(np.sum(num_particles))
        nums = np.insert(np.cumsum(num_particles), 0, 0)
        id = [ ids[nums[i]:nums[i+1]] for i in range(len(num_particles)) ]
            
    header = f.create_group('Header')
    
    header.create_dataset('NumPartTypes', data=np.array(len(num_particles)))
    header.create_dataset('NumParts', data=num_particles)
    
    particles = f.create_group('Particles')
    for i in range(len(num_particles)):
        data = particles.create_group('PartType%d' % i)
        data.create_dataset('Pos', data=pos[i], dtype=np.float64)
        data.create_dataset('Vel', data=vel[i], dtype=np.float64)
        data.create_dataset('Mass', data=mass[i], dtype=np.float64)
        data.create_dataset('ID', data=id[i], dtype=np.int64)

    if time is not None:
        header.create_dataset('Time', data=np.array(time))

    f.close()

