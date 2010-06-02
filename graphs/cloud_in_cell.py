from numpy import uint64, clip, searchsorted, bincount, maximum, minimum, diff, float64, zeros

def cloud_in_cell(x,y, pts_x, pts_y, min_val=False):
    """
    creates a grid of cells with borders x[:] and y[:]
    counts the fraction of particles in each cell, using nearest cell approximation

    x - x borders
    y - y borders
    pts_x - x values of the points
    pts_y - y values of the points
    min_val=False - if True then set the minimum number in each cell as 1, (good for log plots)

    """
    res = zeros((x.size* y.size), dtype=float64)

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


    idx00 = idx_y0 * x.size + idx_x0
    idx01 = idx_y0 * x.size + idx_x1
    idx10 = idx_y1 * x.size + idx_x0
    idx11 = idx_y1 * x.size + idx_x1

    w00 = (1-frac_y) * (1-frac_x)
    w01 = (1-frac_y) * frac_x
    w10 = frac_y * (1-frac_x)
    w11 = frac_y * frac_x

    for idx, weights in zip((idx00,idx01,idx10,idx11), (w00,w01,w10,w11)):
        b = bincount(idx, weights)
        res[:b.size] += b

    if min_val:
        res = maximum(res, 1.0)

    res = res.reshape((y.size, x.size)) / res.sum()

    return res
