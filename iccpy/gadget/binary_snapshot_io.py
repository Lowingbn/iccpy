import numpy as np
from numpy import uint32, uint64, float64, float32

header_names = ('npart', 'mass', 'time', 'redshift', 'flag_sfr', 'flag_feedback', 'npartTotal', 'flag_cooling', \
                'num_files', 'boxsize', 'omega0', 'omegaLambda', 'hubble0', 'flag_stellarage', 'flag_metals', \
                'npartTotalHighWord', 'flag_entropy_instead_u', 'flag_doubleprecision', 'buffer')

header_sizes = ((uint32, 6), (float64, 6), (float64, 1), (float64, 1), (uint32, 1), (uint32, 1), \
             (uint32, 6), (uint32, 1), (uint32, 1), (float64, 1), (float64, 1), (float64, 1), \
             (float64, 1), (uint32, 1), (uint32, 1), (uint32, 6), (uint32, 1), (uint32, 1), (np.uint8, 56))

header_sizes_swap = (('>u4', 6), ('>f8', 6), ('>f8', 1), ('>f8', 1), ('>u4', 1), ('>u4', 1), \
             ('>u4', 6), ('>u4', 1), ('>u4', 1), ('>f8', 1), ('>f8', 1), ('>f8', 1), \
             ('>f8', 1), ('>u4', 1), ('>u4', 1), ('>u4', 6), ('>u4', 1), ('>u4', 1), (np.uint8, 56))

def read_snapshot_file(filename, gas=False, ics=False, cooling=False, accel=False):
    """ Reads a binary gadget file """
    f = open(filename, mode='rb')

    res = {} # dictionary to hold data

    # read the header
    header = read_header(f)

    nparts = header['npart']
    masses = header['mass']

    print 'Particles', nparts
    print 'Masses', masses

    total = nparts.sum()
    print 'Total particles', total
    
    if header['flag_doubleprecision']:
        precision = float64
        print 'Precision: Double'
    else:
        precision = float32
        print 'Precision: Float'

    pos = readu(f, precision, total * 3, header['swap_endian']).reshape((total, 3))
    vel = readu(f, precision, total * 3, header['swap_endian']).reshape((total, 3))

    id = readIDs(f, total, header['swap_endian'])
    
    pmass = []
    #Any particle types which do not have their mass specified in the header will be 
    #found in a mass block
    mass_len = sum([ num for mass, num in zip(masses, nparts) if num>0 and mass==0 ])
    
    if mass_len>0:
        mass_block = readu(f, precision, mass_len, header['swap_endian'])
    
    offset = 0
    for mass, num in zip(masses, nparts):
        if num > 0:
            if mass == 0:
                pmass.append(mass_block[offset:offset+num])
                offset += num
            else:
                pmass.append(mass)
        else:
            pmass.append(0)
    
    if gas:
        ngas = nparts[0]
        res['therm'] = readu(f, float32, ngas, header['swap_endian'])
        if not ics:
            res['rho'] = readu(f, float32, ngas, header['swap_endian'])
            if cooling:
                res['Ne'] = readu(f, float32, ngas, header['swap_endian'])

        if cooling:
            res['NHI'] = readu(f, float32, ngas, header['swap_endian'])
            res['NHeI'] = readu(f, float32, ngas, header['swap_endian'])
            res['NHeIII'] = readu(f, float32, ngas, header['swap_endian'])
                
        if not ics:
            res['sml'] = readu(f, float32, ngas, header['swap_endian'])
        else:
            res['sml'] = 0.0
            
    if accel:
        res['accel'] = readu(f, precision, total * 3, header['swap_endian']).reshape((total, 3))

    f.close()

    res['pos'] = pos
    res['vel'] = vel
    res['mass'] = pmass
    res['id'] = id
    
    return header, res

def write_snapshot_file(filename, header, pos, vel, ids, masses=None, extra_data=None):
    """ Write a binary gadget file (unformatted fortran) """

    # do some checks on the header
    for name, size in zip(header_names, header_sizes):
        if name not in header:
            print 'Missing %s in header file' % name
            raise Exception('Missing %s in header file' % name)

        if np.array(header[name]).size != size[1]:
            msg = 'Header %s should contain %d elements, %d found' % \
                    (name, size[1], np.array(header[name]).size) 
            print msg
            raise Exception(msg)

    # do some checks on the data
    nparts = header['npart']
    
    mass_len = sum([ num for num, mass in zip(nparts, header['mass']) if mass==0 ])
    if mass_len!=0 and mass_len!=len(masses):
        raise Exception('bad mass values')

    #should we write in single or double precision
    if header['flag_doubleprecision']==0:
        precision = float32
    else:
        precision = float64

    # ok so far, lets write the file
    f = open(filename, 'wb')

    # write the header
    np.array(256, uint32).tofile(f)

    for name, size in zip(header_names, header_sizes):
        np.array(header[name]).astype(size[0]).tofile(f)
    # final block
    np.array(256, uint32).tofile(f)
    
    # write the data
    writeu(f, pos.astype(precision))
    writeu(f, vel.astype(precision))
    writeu(f, ids)
    
    if mass_len>0:
        writeu(f, masses.astype(precision))
        
    if extra_data is not None:
      for data in extra_data:
          writeu(f, data.astype(precision))
 
    # all done!
    f.close()

def read_snapshot_header(filename):
    """ Reads just the header from a binary GADGET file """
    f = open(filename, mode='rb')
    header = read_header(f)
    f.close()
    
    return header

def read_header(f):
    """ Read the binary GADGET header file into a dictionary """
    block_size = np.fromfile(f, uint32, 1)[0] 
    if  block_size == 256: # 256 byte header
        header = dict(((name, np.fromfile(f, dtype=size[0], count=size[1])) \
                       for name, size in zip(header_names, header_sizes)))
        header['swap_endian'] = False

        assert(np.fromfile(f, uint32, 1)[0] == 256)
    else:
        header = dict(((name, np.fromfile(f, dtype=size[0], count=size[1])) \
                       for name, size in zip(header_names, header_sizes_swap)))
        header['swap_endian'] = True

        assert(np.fromfile(f, '>u4', 1)[0] == 256)
    return header
    
def readIDs(f, count=None, swap_endian=False):
    """ Read a the ID block from a binary GADGET snapshot file """
    data_size = np.fromfile(f, rtype(uint32, swap_endian), 1)[0]
    
    count = int(count)        
    if data_size / 4 == count: dtype = uint32
    elif data_size / 8 == count: dtype = uint64
    else: raise Exception('Incorrect number of IDs requested')
    
    print "ID size: ", dtype

    ask_size = np.dtype(dtype).itemsize * count
    if ask_size > data_size:
        raise Exception('Data requested larger than buffer')

    arr = np.fromfile(f, rtype(dtype, swap_endian), count)
    final_block = np.fromfile(f, rtype(uint32, swap_endian), 1)[0]
    
    print final_block

    # check the flag at the beginning corresponds to that at the end
    assert(data_size == final_block)

    return arr   

def readu(f, dtype=None, count=None, swap_endian=False):
    """ Read a numpy array from the unformatted fortran file f """  
    data_size = np.fromfile(f, rtype(uint32, swap_endian), 1)[0]
    read_size = data_size

    if dtype is not None:
        count = int(count)    
        ask_size = np.dtype(dtype).itemsize * count
        if ask_size > read_size:
            raise Exception('Data requested larger than buffer')

        arr = np.fromfile(f, rtype(dtype, swap_endian), count)
        read_size = read_size - ask_size
    else:
        arr = None

    f.seek(read_size, 1)
    final_block = np.fromfile(f, rtype(uint32, swap_endian), 1)[0]

    # check the flag at the beginning corresponds to that at the end
    assert(data_size == final_block)

    return arr

def writeu(f, arr=None):
    """ Write a numpy array to the unformatted fortran file f """
    if arr is None:
        np.array(0, dtype=np.uint64).tofile(f)
        return

    data_size = arr.size * arr.dtype.itemsize
    data_size = np.array(data_size, dtype=uint32)

    data_size.tofile(f)

    if arr is not None:
        arr.tofile(f)

    data_size.tofile(f)

def rtype(t, swap_endian):
    if not swap_endian: 
        return t
    else:
        dt = np.dtype(t)
        return dt.newbyteorder('>')
