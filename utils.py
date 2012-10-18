from numpy import *

class Struct:
    """
    Generic class that converts named constructor arguments into class members
    """
    def __init__(self, *args, **kwds):
        for arg in args:
            self.__dict__.update(arg)
        self.__dict__.update(kwds)
        
    def __repr__(self):
        return str(self.__dict__)

def match(arr1, arr2, arr2_sorted=False):
    """
    For each element in arr1 return the index of the element with the
    same value in arr2, or -1 if there is no element with the same value.
    Setting arr2_sorted=True will save some time if arr2 is already sorted
    into ascending order.

    It is assumed that each element in arr1 only occurs once in arr2.
    """

    # Workaround for a numpy bug (<=1.4): ensure arrays are native endian
    # because searchsorted ignores endian flag
    if not(arr1.dtype.isnative):
        arr1_n = asarray(arr1, dtype=arr1.dtype.newbyteorder("="))
    else:
        arr1_n = arr1
    if not(arr2.dtype.isnative):
        arr2_n = asarray(arr2, dtype=arr2.dtype.newbyteorder("="))
    else:
        arr2_n = arr2

    # Sort arr2 into ascending order if necessary
    tmp1 = arr1_n
    if arr2_sorted:
        tmp2 = arr2_n
        idx = slice(0,len(arr2_n))
    else:
        idx = argsort(arr2_n)
        tmp2 = arr2_n[idx]

    # Find where elements of arr1 are in arr2
    ptr  = searchsorted(tmp2, tmp1)

    # Make sure all elements in ptr are valid indexes into tmp2
    # (any out of range entries won't match so they'll get set to -1
    # in the next bit)
    ptr[ptr>=len(tmp2)] = 0
    ptr[ptr<0]          = 0

    # Return -1 where no match is found
    ind  = tmp2[ptr] != tmp1
    ptr[ind] = -1

    # Put ptr back into original order
    ind = arange(len(arr2_n))[idx]
    ptr = where(ptr>= 0, ind[ptr], -1)

    return ptr
    
def linear_index_wts(grid_pts, pts):
    """ 
    find the two grit_pts indices for each point and a weight such that the linear interpolation can be found by
    f(pts) = w * grid_pts[idx0] + (1-w) * grid_pts[idx1]

    returns w, idx0, idx1

    Note, outside the boundaries we set to the edge value (i.e w=0 or 1)
    """
    idx = searchsorted(grid_pts, pts)
    idx0 = maximum(idx-1, 0)
    idx1 = minimum(idx, grid_pts.size-1)
    
    dx = diff(grid_pts)
    w = 1.0 - (pts - grid_pts[idx0]) / dx[minimum(idx0, grid_pts.size-2)]
    
    w = where(equal(idx0, idx1), 1.0, w)

    return w, idx0, idx1

def interp3d(x_grid, y_grid, z_grid, vals, x_pts, y_pts, z_pts):
    """ 3d linear interpolation """
    wx, x_idx0, x_idx1 =  linear_index_wts(x_grid, x_pts)
    wy, y_idx0, y_idx1 =  linear_index_wts(y_grid, y_pts)
    wz, z_idx0, z_idx1 =  linear_index_wts(z_grid, z_pts)

    # interpolate 8 grid points
    res = wx * (wy * (wz * vals[x_idx0, y_idx0, z_idx0] + (1-wz) * vals[x_idx0, y_idx0, z_idx1]) +
        (1-wy) * (wz * vals[x_idx0, y_idx1, z_idx0] +  (1-wz) * vals[x_idx0, y_idx1, z_idx1])) +\
        (1-wx) * (wy * (wz * vals[x_idx1, y_idx0, z_idx0] + (1-wz) * vals[x_idx1, y_idx0, z_idx1]) +
        (1-wy) * (wz * vals[x_idx1, y_idx1, z_idx0] +  (1-wz) * vals[x_idx1, y_idx1, z_idx1]))

    return res  
    
import numpy as np

def points_uniform_on_sphere(n):
    """ Generate a uniform distribution of n points on the surface of a unit sphere """
    pts = np.empty([n, 3])

    inc = np.pi * (3 - np.sqrt(5))
    off = 2.0 / n
    for k in range(0, n):
        y = k * off - 1 + (off / 2)
        r = np.sqrt(1 - y*y)
        phi = k * inc
        pts[k] = np.array([np.cos(phi)*r, y, np.sin(phi)*r])

    return pts
