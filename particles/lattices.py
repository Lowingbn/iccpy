from numpy import mgrid, empty, float64, sqrt, array, pi


def ccp_lattice(nx,ny,nz):
  """
  creates a cubic close packed (ccp) lattice

  A CCP lattice has points at the corners and the centres of faces of a cube
  Using a unit grid the nearest neighbours are sqrt(0.5) apart

  returns positions (3,N) array
  """
  pos_base = mgrid[:nx,:ny,:nz]

  # vectors for offsets

  offset1 = array((0.5, 0.5, 0.0), dtype=float64)
  offset2 = array((0.5, 0.0, 0.5), dtype=float64)
  offset3 = array((0.0, 0.5, 0.5), dtype=float64)


  # create the 4 layers
  pos = empty((4,3,nx,ny,nz), dtype=float64)
  pos[0] = pos_base
  pos[1] = pos[0] + offset1.reshape((3,1,1,1))
  pos[2] = pos[0] + offset2.reshape((3,1,1,1))
  pos[3] = pos[0] + offset3.reshape((3,1,1,1))

  pos = pos.swapaxes(1,-1)
  pos = pos.reshape((4*nx*ny*nz,3))
  return pos
  
def hcp_lattice(nx,ny,nz):
  """ 
  creates a hexagonal close packed (hcp) lattice
  """
  pos_base = mgrid[:nx,:ny,:nz]

  # vectors for multiplies and offsets
  scales = array((1, sqrt(3), (2.0/3.0)* sqrt(6)), dtype=float64)
  offset1 = array((0.5, sqrt(3)*0.5, 0.0), dtype=float64)
  offset2 = array((0.5, sqrt(3)/6.0, sqrt(6.0)/3.0), dtype=float64)

  # create the 4 layers
  pos = empty((4,3,nx,ny,nz), dtype=float64)
  pos[0] = pos_base * scales.reshape((3,1,1,1))
  pos[1] = pos[0] + offset1.reshape((3,1,1,1))
  pos[2:] = pos[:2] + offset2.reshape((1,3,1,1,1))

  pos = pos.swapaxes(1,-1)
  pos = pos.reshape((4*nx*ny*nz,3))
  return pos, scales * (nx,ny,nz)


if __name__=='__main__':

  pos, box_size = hcp_lattice(10,10,10)
  print box_size
  print 'number of points', pos.shape[0]
  vol = box_size[0]*box_size[1] * box_size[2]
  print 'volume', vol
  print 'density', pos.shape[0] / vol
  print 'packing factor', pos.shape[0] * (pi / 6.0) / vol

  pos = ccp_lattice(5,10,15)
  
  import pylab as pl
  pl.plot(pos[:,0], pos[:,2], 'bx')
  pl.show()
  

