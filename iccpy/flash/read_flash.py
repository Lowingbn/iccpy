from tables import openFile
from numpy import mgrid, add, multiply, ones, flatnonzero, empty, float64, repeat

def read_flash_names(filename, names):
    """ read arrays of the given names in a flash file """
    f = openFile(filename, 'r')

    res = {}

    params = f.getNode(f.root, 'simulation parameters')
    res['time'] = params[0][1]

    for name in names:
        res[name] =  f.getNode(f.root, name).read()
    f.close()
    return res

def read_cell_centres(filename):
    """ get cell centres of flash file cells """
    names = ['block size', 'coordinates']
    res = read_flash_names(filename, names) 
  
    block_size =  res['block size']
    coords = res['coordinates']

    if coords.shape[1]==3:
        z,y,x = (mgrid[:8,:8,:8]-3.5) * (1.0/8.0)
 
        x = multiply.outer(block_size[:,0], x) + multiply.outer(coords[:,0], ones(x.shape, dtype=block_size.dtype))
        y = multiply.outer(block_size[:,1], y) + multiply.outer(coords[:,1], ones(y.shape, dtype=block_size.dtype))
        z = multiply.outer(block_size[:,2], z) + multiply.outer(coords[:,2], ones(z.shape, dtype=block_size.dtype))
  
        return x,y,z

    if coords.shape[1]==2:
        y,x = (mgrid[:8,:8]-3.5) * (1.0/8.0)
 
        x = multiply.outer(block_size[:,0], x) + multiply.outer(coords[:,0], ones(x.shape, dtype=block_size.dtype))
        y = multiply.outer(block_size[:,1], y) + multiply.outer(coords[:,1], ones(y.shape, dtype=block_size.dtype))

        return x,y
    x = (mgrid[:8]-3.5) * (1.0/8.0)
    x = multiply.outer(block_size[:,0], x) + multiply.outer(coords[:,0], ones(x.shape, dtype=block_size.dtype))
    
    return (x,)

def read_runtime_params(filename):
    """ read the cfl value """
    f = openFile(filename, 'r')
    params = f.getNode(f.root, 'real runtime parameters').read()

    f.close()

    params = dict((name.rstrip(), val) for name, val in params)

    return params

def get_data(filename, vars=('dens', 'pres', 'lamb', 'temp', 'eint')):
    """ 
    get_data(filename, vars=('dens', 'pres', 'lamb', 'temp', 'eint')):
    read flash data for a given filename
    """
    coords = read_cell_centres(filename)
    dim = len(coords)
    names = ['node type', 'block size'] + list(vars)
    data = read_flash_names(filename, names) 

    leaf_blocks = flatnonzero(data['node type']==1)

    result = {}
    for var in vars:
        result[var] = data[var][leaf_blocks].flatten()

    result['cell size'] = repeat(data['block size'][leaf_blocks, 0] / 8.0, 8**dim)

    pos = empty((leaf_blocks.size*(8**dim), dim), dtype=float64)

    for i in range(dim):
        pos[:,i] = coords[i][leaf_blocks].flatten()


    result['pos'] = pos
    result['time'] = data['time']

    return result

if __name__=='__main__':
    filename = '../Halo_hdf5_chk_0000'
    x,y,z = read_cell_centres(filename)
    names = ['node type', 'dens']
    res = read_flash_names(filename, names) 

    dens = res['dens']
    node_type = res['node type']

    # keep only the leaf blocks
    indices = flatnonzero(node_type==1)
  
    dens = dens[indices].flatten()
    x    =  x[indices].flatten()
    y    =  y[indices].flatten()
    z    =  z[indices].flatten()

    # take a column of data
    xmin, xmax = 3.0e25, 3.1e25
    indices = flatnonzero((z>xmin) & (z<xmax) & (x>xmin) & (x<xmax))
 
    dens = dens[indices]
    x    =  x[indices]
    y    =  y[indices]
    z    =  z[indices]

    from pylab import plot, show, log
    plot(y, log(dens), 'rx')
    show()
    #print unique(z)


