from tables import openFile
from pylab import plot, show, loglog, xlabel, ylabel, legend
from read_flash import read_flash_names, read_cell_centres
import numpy as np

yr = 3.155e7 # year in seconds
kpc = 3.0857e21 # kpc in cm 
G = 6.674e-8
m_p = 1.6726e-24 # proton mass in g

def dens_profile(filename, centre, name='dens'):
  x,y,z = read_cell_centres(filename)
  names = ['node type', name, 'block size', 'pres']
  res = read_flash_names(filename, names) 

  var = res[name]
  pres = res['pres']
  node_type = res['node type']
  time = res['time']
  bs = res['block size']

  from numpy import flatnonzero, unique, square, sqrt
  # keep leaf nodes
  indices = flatnonzero(node_type==1)


  var  = var[indices].flatten()
  pres = pres[indices].flatten()
  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()
  cell_size   =  bs[indices,0]/8.0/kpc

#  if name=='dens':
#    print var.shape
#    print cell_size.shape
#    print 'Total mass %e'% float((cell_size**3 * var).sum())
  print 'cell size in', cell_size.min(), cell_size.max(),'^3 kpc'
  
  
  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]

  rad = sqrt(square(x) + square(y) + square(z))
  
  return rad, var, time

def mass_profile(filename, centre):
  x,y,z = read_cell_centres(filename)
  names = ['node type', 'dens', 'block size']
  res = read_flash_names(filename, names) 

  dens = res['dens']
  node_type = res['node type']
  time = res['time']
  bs = res['block size']

  from numpy import flatnonzero, unique, square, sqrt
  # keep leaf nodes
  indices = flatnonzero(node_type==1)


  dens  = dens[indices].flatten()
  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()
  cell_size   =  bs[indices,0]/8.0
  cell_size = np.repeat(cell_size, 8**3)

  mass = ((dens * cell_size) * cell_size) * cell_size

  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]

  rad = sqrt(square(x) + square(y) + square(z))

  idx = np.argsort(rad)
  rad = rad[idx]
  mass = np.cumsum(mass[idx])  
  
  spots = np.flatnonzero((np.diff(rad)/ rad[1:])>0.01)
  rad = rad[spots]
  mass = mass[spots]

  return rad, mass, time


def cool_profile(filename, centre):
  x,y,z = read_cell_centres(filename)
  names = ['node type', 'dens', 'block size', 'lamb', 'H   ', 'eint']
  res = read_flash_names(filename, names) 

  dens = res['dens']
  h = res['H   ']
  eint = res['eint'] # specific internal energy
  lamb = res['lamb']
  node_type = res['node type']
  time = res['time']

  from numpy import flatnonzero, unique, square, sqrt
  # keep leaf nodes
  indices = flatnonzero(node_type==1)
  
  
  dens = dens[indices].flatten()
  eint = eint[indices].flatten()
  lamb = lamb[indices].flatten()
  h    = h[indices].flatten()
  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()

  
  n_h = (dens * h) / m_p
  cool = -lamb * n_h * n_h
  t_cool = (eint * dens)/ cool 
  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]

  rad = sqrt(square(x) + square(y) + square(z))
  print lamb.max()
  return rad, t_cool, time

def vel_profile(filename, centre):
  """ get the velocity profile of the gas """
  x,y,z = read_cell_centres(filename)
  names = ['node type', 'velx', 'vely', 'velz']
  res = read_flash_names(filename, names) 

  velx = res['velx']
  vely = res['vely']
  velz = res['velz']
  
  node_type = res['node type']
  time = res['time']

  from numpy import flatnonzero, unique, square, sqrt
  # keep leaf nodes
  indices = flatnonzero(node_type==1)
  
  
  velx = velx[indices].flatten()
  vely = vely[indices].flatten()
  velz = velz[indices].flatten()
  
  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()

  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]

  rad = sqrt(square(x) + square(y) + square(z))

  # radial vector
  rx, ry, rz = x / rad, y/rad, z/rad
  # radial component of velocity
  vel = velx * rx + vely * ry + velz * rz  
  
  return rad, vel, time

def entropy_profile(filename, centre):
  """ radial entropy profile """
  x,y,z = read_cell_centres(filename)
  names = ['node type', 'dens', 'block size', 'pres']

  res = read_flash_names(filename, names) 

  dens = res['dens']
  node_type = res['node type']
  time = res['time']
  pres = res['pres']

  from numpy import flatnonzero, unique, square, sqrt
  # keep leaf nodes
  indices = flatnonzero(node_type==1)
  
  
  dens = dens[indices].flatten()
  pres = pres[indices].flatten()
  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()

  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]


  rad = sqrt(square(x) + square(y) + square(z))
  
  entropy = pres * (dens**-1.6666667)
  return rad, entropy, time

def gravity_profile(filename, centre, n=100):
  """ get the gravity profile of the gas, n=num of rad pts"""

  rad, gpot, time = dens_profile(filename, centre, 'gpot')


  idx = np.argsort(rad)
  rad = rad[idx]
  gpot = gpot[idx]

  log_rg = np.log(rad.max()/rad.min())
  keep = (np.log(rad) - np.log(rad.min())) * (n / log_rg)
  keep = keep.astype(np.int32)
  print keep.shape
  idx = np.flatnonzero(np.diff(keep))
  print idx.shape
  
  #idx = np.array([0] + list(idx), dtype=np.int32)
  rad = rad[idx]

  av_pot = np.empty_like(rad) 
  av_pot[1:] = np.diff(np.cumsum(gpot)[idx]) / np.diff(idx)
  av_pot[0] = gpot[idx[0]]

  
  gpot = av_pot
  
  grav = np.zeros_like(gpot)
  grav[1:] = np.diff(gpot)/np.diff(rad)

  return rad, grav, time

def gradp_profile(filename, centre):
  """ get the pressure gradient of the gas """

  rad, pres, time = dens_profile(filename, centre, 'pres')

  n = 500 # number of graph points
  import numpy as np
  idx = np.argsort(rad)
  rad = rad[idx]
  pres = pres[idx]

  log_rg = np.log(rad.max()/rad.min())
  keep = (np.log(rad) - np.log(rad.min())) * (n / log_rg)
  keep = keep.astype(np.int32)
  print keep.shape
  idx = np.flatnonzero(np.diff(keep))
  print idx.shape
  
  rad = rad[idx]

  av_pres = np.empty_like(rad) 
  av_pres[1:] = np.diff(np.cumsum(pres)[idx]) / np.diff(idx)
  av_pres[0] = pres[idx[0]]

  
  pres = av_pres
  
  gradp = np.zeros_like(pres)
  gradp[1:] = np.diff(pres)/np.diff(rad)

  return rad, gradp, time

def rho_grav_profile(filename, centre):
  """ get the rho*grad Phiprofile """

  rad, grav, time = gravity_profile(filename, centre)

  rad2, dens, time = dens_profile(filename, centre)

  import numpy as np
  idx = np.searchsorted(rad, rad2)
  idx = np.clip(idx, 0, grav.size-1)
  grav = grav[idx]
  rho_grav = dens * grav

 
  return rad2, rho_grav, time

def cold_gas(filenames, eint_cutoff):
  """ find the amount of cold gas in each file """
  times = np.zeros(len(filenames), dtype=np.float64)
  masses = np.zeros_like(times)
  from numpy import flatnonzero, unique, square, sqrt

  for i, filename in enumerate(filenames):
    names = ['node type', 'dens', 'block size', 'eint']
    res = read_flash_names(filename, names) 

    dens = res['dens']
    eint = res['eint']
    node_type = res['node type']
    times[i] = res['time']
    block_size = res['block size']


    # keep leaf nodes
    indices = flatnonzero(node_type==1)

    dens = dens[indices].flatten()
    eint = eint[indices].flatten()

    width = block_size[indices,0] / 8.0 # width in cm
    width = np.repeat(width, 8**3)
    mass = (dens*width) * width * width

    indices = np.flatnonzero(eint < eint_cutoff)
    if indices.size > 0:
      masses[i] = mass[indices].sum()

  times /= yr * 1e9 # convert to Gyr
  return times, masses

def cell_size_profile(filename, centre):
  x,y,z = read_cell_centres(filename)
  names = ['node type', 'block size']
  res = read_flash_names(filename, names) 

  node_type = res['node type']
  time = res['time']
  bs = res['block size']

  from numpy import flatnonzero, unique, square, sqrt, repeat
  # keep leaf nodes
  indices = flatnonzero(node_type==1)

  x    =  x[indices].flatten()
  y    =  y[indices].flatten()
  z    =  z[indices].flatten()
  cell_size   =  bs[indices,0]/8.0

  cell_size = repeat(cell_size, 8**3)

  x = x-centre[0]
  y = y-centre[1]
  z = z-centre[2]

  rad = sqrt(square(x) + square(y) + square(z))
  
  return rad, cell_size, time
