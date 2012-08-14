import numpy as np
import h5py
import sys
import iccpy.gadget.binary_snapshot_io

def write_snapshot_file(filename, num_particles, pos, vel, mass, id=None, time=None):
    f = h5py.File('%s.hdf5' % filename, 'w')
    
    if type(num_particles)==int:
        num_particles = np.array([num_particles])
        pos = [ pos ]
        vel = [ vel ]
        mass = [ mass ]
        if id is not None: id = [ id ]
    else:
        num_particles = np.array(num_particles)
    
    if id is None:
        ids = np.arange(np.sum(num_particles), dtype=np.uint64)
        nums = np.insert(np.cumsum(num_particles), 0, 0)
        id = [ ids[nums[i]:nums[i+1]] for i in range(len(num_particles)) ]
            
    header = f.create_group('Header')
    
    header.create_dataset('NumPartTypes', data=np.array([len(num_particles)]),  dtype=np.int32)
    header.create_dataset('NumParts', data=num_particles,  dtype=np.int32)
    
    particles = f.create_group('Particles')
    for i in range(len(num_particles)):
        data = particles.create_group('PartType%d' % i)
        data.create_dataset('Pos', data=pos[i], dtype=np.float64)
        data.create_dataset('Vel', data=vel[i], dtype=np.float64)
        data.create_dataset('Mass', data=mass[i], dtype=np.float64)
        data.create_dataset('ID', data=id[i], dtype=np.uint64)

    if time is not None:
        header.create_dataset('Time', data=np.array(time))

    f.close()

def convert_to_gadget(filename):
    f = h5py.File('%s.hdf5' % filename, 'r')
    
    num_particle_types = f['/Header/NumPartTypes'][0]
    
    if num_particle_types>5:
        print "Error - Gadget only supports 5 dark matter particle types. %s has %d" % (filename, num_particle_types)
    
    num_particles = np.zeros(6)
    pos = []
    vel = []
    mass_block = []
    mass_header = np.zeros(6)
    ids = []

    time = 0
    if 'Time' in f['/Header']: 
        time = f['/Header/Time'][0]
    
    for i in range(num_particle_types):
        group = f['/Particles/PartType%d' % i]
        
        pos.append(np.array(group['Pos'], dtype=np.float64))
        vel.append(np.array(group['Vel'], dtype=np.float64))
        mass = np.array(group['Mass'], dtype=np.float64)
        id = np.array(group['ID'], dtype=np.int64)
        ids.append(id)

        num_particles[i+1] = len(id)

        if len(np.unique(mass))>1 or mass[0]==0:
            #Different particle masses (or zero mass particles) so need to place in mass block
            mass_block.append(mass)
        else:
            mass_header[i+1] = np.unique(mass)[0]
            
    pos = np.concatenate(pos)
    vel = np.concatenate(vel)
    ids = np.concatenate(ids)
                
    if len(mass_block)==0: 
        mass_block = None 
    else:
        mass_block = np.concatenate(mass_block)
    
    header = dict((('num_particles', num_particles),
                   ('mass', mass_header), ('time',time), ('redshift',0), 
                   ('flag_sfr',0) , ('flag_feedback',0), 
                   ('num_particles_total', num_particles), 
                   ('flag_cooling',0), ('num_files',1), ('boxsize',1.0), \
                   ('omega0', 0.3), ('omegaLambda', 0.7), ('hubble0', 1.0), ('flag_stellarage',0), \
                   ('buffer', [0]*56), ('flag_metals', 0), ('npartTotalHighWord', [0,0,0,0,0,0]), 
                   ('flag_entropy_instead_u', 0), ('flag_doubleprecision', 1)))
                
    print "Writing gadget snapshot %s.snp" % filename
    print header
    
    iccpy.gadget.binary_snapshot_io.write_snapshot_file("%s.snp" % filename, header, pos, vel, ids, masses=mass_block)
    
if __name__=="__main__":
    if len(sys.argv)!= 2:
        print "Usage: ./convert.py filename"
        sys.exit()
    
    convert_to_gadget(sys.argv[1])    
