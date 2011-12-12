"""
Code created by Peter Creasey
"""
from tables import openFile
from numpy import mgrid, flatnonzero, empty, float64, repeat, int32, unique, arange, cumprod, array
from itertools import izip

def read_flash_names(f, names):
    """ read arrays of the given names in a flash file """
    res = {}

    # need to check if flash 3
    name = 'simulation parameters'
    if name not in f.root:
        name = 'real scalars'

    params = f.getNode(f.root, name)
    res['time'] = params[0][1]

    for name in names:
        res[name] =  f.getNode(f.root, name).read()

    return res

def blocks_to_cells(block_size, coords, b_shape):
    """
    Find cell coordinates from block coordinates
    block_size - (n,dim) array of block sizes (e.g. bx,by,bz for n blocks)
    coords - (n,dim) array of block centres (e.g cx,cy,cz for n blocks)
    b_shape - shape of the block in cells, e.g. [8,8,8]

    returns cell_pos, cell_size
    (n*M, d) arrays, where M = Prod(b_shape)
    """
    
    dim = block_size.shape[1]
    nb = b_shape

    cells_per_block = cumprod(nb)[-1]

    nblocks = len(block_size)

    if dim==3:
        rev_nb = nb[::-1].reshape(dim,1,1,1)
        p2 = (mgrid[:nb[2],:nb[1],:nb[0]]+0.5*(1.0 - rev_nb)) * (1.0/ rev_nb)
        cell_pos = (p2.reshape((dim,cells_per_block))[::-1]).swapaxes(0,1)
    elif dim==2:
        rev_nb = nb[::-1].reshape(dim,1,1)
        p2 = (mgrid[:nb[1],:nb[0]]+0.5*(1.0 - rev_nb)) * (1.0/ rev_nb)
        cell_pos = (p2.reshape((dim,cells_per_block))[::-1]).swapaxes(0,1)
    else:
        p2 = (arange(nb[0])+0.5*(1.0 - nb[0])) * (1.0/ nb[0])
        cell_pos = (p2.reshape((dim,cells_per_block))[::-1]).swapaxes(0,1)

    all_cell_pos = block_size.reshape((nblocks,1,dim)) * cell_pos + coords.reshape((nblocks,1,dim))
    all_cell_pos = all_cell_pos.reshape((nblocks * cells_per_block, dim))
    all_cell_size = repeat(block_size[:,0] / nb[0], cells_per_block)
    return all_cell_pos, all_cell_size


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
    f = openFile(filename, 'r')
    
    dim,nb = block_shape(f)
    cells_per_block = cumprod(nb)[-1]
    names = ['node type', 'block size', 'coordinates'] + list(vars)
    data = read_flash_names(f, names) 

    
    leaf_blocks = flatnonzero(data['node type']==1)

    result = {}
    for var in vars:
        result[var] = data[var][leaf_blocks].flatten()

    block_sizes = data['block size'][leaf_blocks]
    dim = block_sizes.shape[1]
    
    cell_pos, cell_size = blocks_to_cells(block_sizes, data['coordinates'][leaf_blocks], nb[:dim]) 
    f.close()

    result['cell size'] = cell_size
    result['pos'] = cell_pos
    result['time'] = data['time']

    return result

def block_shape(f):
    """
    find the block shape (nxb, nyb, nzb) given the hdf5 file f
    returns
    dimension, (nxb, nyb, nzb)
    """
    if 'integer scalars' in f.root:    
        params =  f.getNode(f.root, 'integer scalars').read()
        p_dict = dict((name.rstrip(), val) for name, val in params)

        dimension = p_dict['dimensionality']
        nb = empty(dimension, dtype=int32)
        for i,par in enumerate(['nxb', 'nyb', 'nzb'][:dimension]):
            nb[i]= p_dict[par]
    else:
        print dir(f.getNode(f.root, 'block size'))
        dimension = 3
        params = f.getNode(f.root, 'simulation parameters')
        nb = empty(dimension, dtype=int32)
        for i in range(dimension):
            nb[i] = params[0][5+i]

    return dimension, nb

def extent(f, dim):
    """
    return a (dim, 2) array of minimum and maximum value of x,y,z
    """
    
    params =  f.getNode(f.root, 'real runtime parameters').read()
    p_dict = dict((name.rstrip(), val) for name, val in params)

    result = empty((dim,2), dtype=float64)

    for i in range(dim):
        dir = 'xyz'[i] 
        result[i,0] = p_dict[dir+'min']
        result[i,1] = p_dict[dir+'max']
        
    return result
    
def uniform_grid(filename, vars=['dens']):
    """ quickly get uniform grids for the given variables """
    f = openFile(filename, 'r')
    

    # shape of cells in block (e.g. (8,8,8))
    dim, nb = block_shape(f)

    block_sizes = f.getNode(f.root, 'block size').read()
    nt = f.getNode(f.root, 'node type').read()
    coords = f.getNode(f.root,'coordinates').read()
        

    # only want the leaf blocks
    leaf_blocks = flatnonzero(nt==1)
    coords = coords[leaf_blocks]

    # check it is really a uniform grid
    for i in range(dim):
        bsizes_i = unique(block_sizes[:,i][leaf_blocks])
        #if len(bsizes_i)!=1:
        #    print 'Unique block sizes (dimension %d)'%i,  bsizes_i
        #    raise Exception('Unequal leaf block sizes, not a uniformly adapted mesh')
        
    # assume uniform grid
    dxyz = block_sizes[leaf_blocks[0]]
    
    # index into the uniform grid
    idx = (coords/dxyz).astype(int32)
    
    # shape of the uniform grid
    grid_n = idx.max(0) + 1

    uni_shape = nb * grid_n
    print 'Uniform grid of shape', uni_shape


    result = []
    for var_name in vars:

        flash_vals = f.getNode(f.root, var_name).read()

        val_arr = empty(uni_shape[::-1], dtype= flash_vals.dtype)
        for idx_i,val_block in izip(idx, flash_vals[leaf_blocks]):
            a = nb * idx_i
            b = nb * (idx_i+1)
            
            val_arr[a[2]:b[2], a[1]:b[1], a[0]:b[0]] = val_block

        del flash_vals
        result.append(val_arr)

    min_max = extent(f, dim)
    
    f.close()

    # make the x,y,z arrays
    xyz = []
    for i in range(dim):
        dmin, dmax = min_max[i]
        n_cells = uni_shape[i]
        
        xyz.append((arange(n_cells)+0.5) * ((dmax-dmin)/float(n_cells)) + dmin)
    
    return tuple(xyz), result

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


