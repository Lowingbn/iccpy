from numpy import uint64, clip, searchsorted, bincount, maximum, minimum, diff, float64, zeros, equal, where, zeros_like

def cloud_in_cell(x,y, pts_x, pts_y, pts_vals=None, min_val=False):
    """
    creates a grid of cells with borders x[:] and y[:]
    counts the fraction of particles in each cell, using cloud-in-cell approximation

    x - x borders
    y - y borders
    pts_x - x values of the points
    pts_y - y values of the points
    pts_vals - values for each point (default=None). If set then the grid contains the average value for the cells

    min_val=False - if True then set the minimum number in each cell as 1, (good for log plots)

    """
    density = zeros((x.size*y.size), dtype=float64)
    if pts_vals is not None:
        values = zeros_like(density)

    # find the surrounding cell borders
    idx_x = searchsorted(x, pts_x)
    idx_x0 = maximum(idx_x-1, 0)
    idx_x1 = minimum(idx_x, x.size-1)

    dx = diff(x)[minimum(idx_x0, x.size-2)]
    frac_x = clip((pts_x - x[idx_x0]) / dx, 0, 1)

    idx_y = searchsorted(y, pts_y)
    idx_y0 = maximum(idx_y-1, 0)
    idx_y1 = minimum(idx_y, y.size-1)

    dy = diff(y)[minimum(idx_y0, y.size-2)]
    frac_y = clip((pts_y - y[idx_y0]) / dy, 0, 1)

    outside_box = equal(idx_x,0) | equal(idx_x,x.size) | equal(idx_y, 0) | equal(idx_y, y.size)
 

    idx00 = idx_y0 * x.size + idx_x0
    idx01 = idx_y0 * x.size + idx_x1
    idx10 = idx_y1 * x.size + idx_x0
    idx11 = idx_y1 * x.size + idx_x1

    w00 = where(outside_box, 0, (1-frac_y) * (1-frac_x))
    w01 = where(outside_box, 0, (1-frac_y) * frac_x)
    w10 = where(outside_box, 0, frac_y * (1-frac_x))
    w11 = where(outside_box, 0, frac_y * frac_x)

    for idx, weights in zip((idx00,idx01,idx10,idx11), (w00,w01,w10,w11)):
        b = bincount(idx, weights)

        density[:b.size] += b

        if pts_vals is not None:
            b = bincount(idx, weights * pts_vals)
            values[:b.size] += b


    res = density

    if min_val:
        res = maximum(res, 1.0)

    if pts_vals is not None:
        res = values / res
        if min_val:
            res = maximum(res, pts_vals.min())
    else:
        # normalise
        res = res / res.sum()


    res = res.reshape((y.size, x.size))

    return res
