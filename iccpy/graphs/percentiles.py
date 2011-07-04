"""
Peter Creasey 2011
p.e.creasey.00@gmail.com

Histogram code with percentiles
"""
from numpy import sort, array, int32, empty, arange,equal, flatnonzero, digitize
        
def percentiles(data, tiles=(0.15, 0.85)):
    """ find the percentiles of the data (default 0.15, 0.85) """
    sdata = sort(data)
    size = len(data) - 1
    pos = array(tiles) * size
    idx0 = pos.astype(int32)
    idx1 = idx0+1
    wts = idx1 -pos
    return wts * sdata[idx0] + (1.0-wts) * sdata[idx1]

    
def p_histogram(x, y, bins=None, bin_range=None, perc=(0.15,0.85)):
    """
    histogram of the percentiles of the data
    returns arrays vals,bins
    vals.shape=(len(bins)-1, len(perc)+1)
    vals[i,0] = mean of bin i
    vals[i,1:] = percentiles of bin i
    
    
    """
    # check if we were given the bins explicitly, or just a number and range
    try:
        b = int(bins)
        if bin_range is None:
            bin_range = (x.min(), x.max())
        bins = arange(b+1) * (bin_range[1]-bin_range[0])/float(b) + bin_range[0]
    except TypeError:
        # bins is a sequence
        bins = array(bins)
        

    # find which data is in which bin
    x_bin = digitize(x, bins)

    # 1 less value than bin edges, store the mean in vals[:,0]
    vals = empty((len(bins)-1, len(perc)+1), dtype=y.dtype)
    
    for i in range(len(bins)-1):
        idx = flatnonzero(equal(x_bin,i+1))

        if idx.size==0:
            vals[i,:] = 0
            continue

        vals[i,1:] = percentiles(y[idx], perc)
        vals[i,0] = y[idx].mean()
            

    return vals, bins
