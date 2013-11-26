from numpy import cumprod, array


def readf(filename, items):
    """ 
    Read IDL arrays from a file 
    
    filename - path of the file
    items - iterable of (func, shape) where func is applied to each split string of the file, 
            then made into an array

    e.g.
    snaps, vels = readf('sfr.dat', [(int, 131), (float, (131,3))])

    """

    fl = open(filename, 'r')
    data = ' '.join(fl.readlines()).split()
    fl.close()
    
    i0 = 0
    for dfunc, shape in items:
        nvals = cumprod(shape)[-1]
        arr = array([dfunc(x) for x in data[i0:i0+nvals]]).reshape(shape)
        yield arr
        i0 += nvals

    return
