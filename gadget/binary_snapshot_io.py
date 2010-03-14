import numpy as np
from numpy import uint32, uint64, float64, float32

mpc = 3.0857e24 # cm 
km = 1.0e5 # cm
msun = 2.0e33 # solar mass in g

header_names = ('npart', 'mass', 'time', 'redshift', 'flag_sfr', 'flag_feedback', 'npartTotal', 'flag_cooling', \
               'num_files', 'boxsize', 'omega0', 'omegaLambda', 'hubble0', 'flag_stellarage', 'flag_metals', \
               'npartTotalHighWord', 'flag_entropy_instead_u', 'flag_doubleprecision', 'buffer')

header_sizes =((uint32,6), (float64, 6), (float64, 1), (float64, 1), (uint32, 1), (uint32,1), \
       (uint32, 6), (uint32,1), (uint32,1), (float64, 1), (float64, 1), (float64, 1), \
       (float64, 1), (uint32,1), (uint32, 1), (uint32, 6), (uint32, 1), (uint32, 1), (np.uint8, 56))


def read_snapshot_file(filename, physical=False, ics=False, cooling=False):
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
  else:
    precision = float32

  pos = readu(f, precision, total*3).reshape((total, 3))
  vel = readu(f, precision, total*3).reshape((total, 3))

  id = readIDs(f, total)
  
  pmass = []
  for mass, num in zip(masses, nparts):
    if num > 0:
      if mass==0:
        pmass.append(readu(f, precision, num))
      else:
        pmass.append(mass)
    else:
      pmass.append(0)

  ngas = nparts[0]
  
  if ngas > 0:
    res['therm'] = readu(f, float32, ngas)
    if not ics:
       res['rho'] = readu(f, float32, ngas)
       if cooling:
         res['Ne'] = readu(f, float32, ngas)

    if cooling:
        res['NHI']  = readu(f, float32, ngas)
        res['NHeI'] = readu(f, float32, ngas)
        res['NHeIII'] = readu(f, float32, ngas)
        
    if not ics:
        res['sml']  =  readu(f, float32, ngas)
    else:
        res['sml'] = 0.0


  h = header['hubble0']
  z = header['redshift']

  lunit = 1.0/h # Mpc
  munit = 1.e10/h
  vunit = np.sqrt(1.0/(1+z))*1.e5/km          # km/s

  print ' units: I assume that lunit=Mpc/h, munit=10^10 modot0/h, v_unit=km/s '

  if physical:
    print ' converting to physical Mpc, km/s, M_odot, g/cm^3'
    pos   = pos * lunit/(1 + z)
    if ngas > 0:
        res['sml']   = res['sml'] * lunit/(1.+z)
        res['rho']   = res['rho'] * munit / (lunit * mpc / (1 + z))**3.0 * msun
        res['therm'] = res['therm'] * vunit**2 * (1 + z)

    vel   = vel * vunit/np.sqrt(1 + z)
    mass  = mass * munit
    pmass = [p * munit for p in pmass]

    munit = msun
    lunit = mpc
    vunit = km

  f.close()


  res['pos'] = pos
  res['vel'] = vel
  res['mass'] = pmass
  res['id'] = id
  
  return header, res

def write_snapshot_file(filename, header, data, physical=True, cooling=False, ics=False, gamma=1.66667):
  """ Write a binary gadget file (unformatted fortran) """

  # do some checks on the header
  for name, size in zip(header_names, header_sizes):
    if name not in header:
      print 'Missing %s in header file'%name
      raise Exception('Missing %s in header file'%name)

    if np.array(header[name]).size != size[1]:
      msg = 'Header %s should contain %d elements, %d found'%\
          (name, size[1], np.array(header[name]).size) 
      print msg
      raise Exception(msg)


  # do some checks on the data
  nparts = header['npart']
  total = np.array(nparts).sum()
  ngas = nparts[0]

  for num, mass, pmass in zip(nparts, header['mass'], data['mass']):
    if num > 0 and mass==0:
      if len(pmass) != num:
        raise Exception('bad mass values')

  init = [('pos', float32, total*3), ('vel', float32, total*3), \
           ('id', uint32, total)]
  a = []
  if ngas > 0:
    a.append(('therm', float32, ngas))
    
    if not ics:
      a.append(('rho', float32, ngas))
      
      if cooling:
        a.append(('Ne', float32, ngas))

    if cooling:
      a.append(('NHI', float32, ngas))
      a.append(('NHeI', float32, ngas))
      a.append(('NHeIII', float32, ngas))

  if not ics:
    a.append(('sml', float32, ngas))


  for name, dtype, size in a:
    if data[name].size != size:
      raise Exception('%s has wrong size'%name)

  # ok so far, lets write the file
  f = open(filename, 'wb')


  if physical:
    # convert to co-moving
    print ' converting to physical Mpc, km/s, M_odot, g/cm^3'
    h = header['hubble0']
    z = header['redshift']

    lunit = 1.0/h # Mpc
    munit = 1.e10/h
    vunit = np.sqrt(1.0/(1+z))*1.e5/km          # km/s

    data['pos']   = data['pos'] * (1+z)/ lunit
    if ngas > 0:
        data['sml']   = data['sml'] / (lunit/(1.+z))
        data['rho']   = data['rho'] / (munit / (lunit * mpc / (1 + z))**3.0 * msun)
        if header['flag_entropy_instead_u']==0:
          data['therm'] = data['therm'] / (vunit**2 * (1 + z))
        else:
          data['therm'] = data['therm'] / (vunit**2 * (1 + z)) * (munit / (lunit * mpc / (1 + z))**3.0 * msun)**(gamma-1)

    data['vel']   = data['vel'] / (vunit/np.sqrt(1 + z))
    header['mass']  = [m / munit for m in header['mass']]
    data['mass'] = [p / munit for p in data['mass']]

  # write the header
  np.array(256, uint32).tofile(f)

  for name, size in zip(header_names, header_sizes):
    np.array(header[name]).astype(size[0]).tofile(f)
  # final block
  np.array(256, uint32).tofile(f)
  
  # write the data
  for name, dtype, size in init:
    writeu(f, data[name].astype(dtype))

  # the masses
  for mass, pmass in zip(header['mass'], data['mass']):
    if mass==0:
      writeu(f, pmass.astype(float32))

  for name, dtype, size in a:
    writeu(f, data[name].astype(dtype))
 
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
  assert(np.fromfile(f, uint32, 1)[0]==256) # 256 byte header


  header = dict(((name, np.fromfile(f, dtype=size[0], count=size[1])) \
                     for name, size in zip(header_names, header_sizes)))


  assert(np.fromfile(f, uint32, 1)[0]==256)
  return header
  
def readIDs(f, count=None):
  """ Read a the ID block from a binary GADGET snapshot file """
  data_size = np.fromfile(f, uint32, 1)[0]
  
  count = int(count)    
  if data_size/4==count: dtype=uint32
  elif data_size/8==count: dtype=uint64
  else: raise Exception('Incorrect number of IDs requested')

  ask_size = np.dtype(dtype).itemsize * count
  if ask_size > data_size:
    raise Exception('Data requested larger than buffer')

  arr = np.fromfile(f, dtype, count)
  final_block = np.fromfile(f, uint32, 1)[0]

  # check the flag at the beginning corresponds to that at the end
  assert(data_size==final_block)

  return arr  

def readu(f, dtype=None, count=None):
  """ Read a numpy array from the unformatted fortran file f """
  data_size = np.fromfile(f, uint32, 1)[0]
  read_size = data_size

  if dtype is not None:
    count = int(count)  
    ask_size = np.dtype(dtype).itemsize * count
    if ask_size > read_size:
      raise Exception('Data requested larger than buffer')

    arr = np.fromfile(f, dtype, count)
    read_size = read_size - ask_size
  else:
    arr=None

  f.seek(read_size, 1)
  final_block = np.fromfile(f, uint32, 1)[0]

  # check the flag at the beginning corresponds to that at the end
  assert(data_size==final_block)

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
