"""
Module for making scatter plots where the points are coloured by the local density of points
Peter Creasey - 2013
"""
from numpy import flatnonzero, empty, float64, log
from scipy.spatial import cKDTree
def density_colours(x, y, xlims, ylims):
    """ 
    find colours for each point using a density estimate, normalise to 0-1     
    Usage: (x,y,xlims, ylims)

    x/ y        - numpy arrays of coordinates
    xlims/ylims - tuples of (min, max), used to cut the particles outside

    returns indices, clrs
    
    indices - the indices of the points within the limits
    clrs    - values in [0,1] correspoding to local density

    Example:
    idx, clr = density_colours(log10(density), log10(temp), [-25,-21], [4,9])
    pl.scatter(density[idx], temp[idx], s=0.25, c=clr, edgecolors='none', cmap='jet')
    """
    assert(len(x)==len(y))
    assert(len(xlims)==2)
    assert(len(ylims)==2)

    # only plot the points within the limits
    idx = flatnonzero((y.ravel() > ylims[0]) &  (y.ravel() < ylims[1]) & 
                      (x.ravel() > xlims[0]) &  (x.ravel() < xlims[1]))

    

    pos = empty((len(idx), 2), dtype=float64)
    pos[:,0] = x.ravel()[idx]*(ylims[1]-ylims[0])
    pos[:,1] = y.ravel()[idx]*(xlims[1]-xlims[0])
    
    # Find the distance to the 10th nearest neighbour as an estimate of local density
    # (10 should be ok as 2d)
    tree = cKDTree(pos)
    nbr_dist = tree.query(pos, 10)[0][:,9]

    c = -log(nbr_dist)


    c = (c-c.min())/(c.max() - c.min())

    return idx, c

        
