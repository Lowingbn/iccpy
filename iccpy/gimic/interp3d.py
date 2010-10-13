from numpy import searchsorted, diff, maximum, minimum, where, equal

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
    
