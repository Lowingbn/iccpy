from numpy import where, pi

def sph_kernel(w, h):
    """ 
    cubic spline kernel as used in  Gadget-2 
    usable with numpy arrays
    """
    r = w/h

    res1 = 2 * (1-r)*(1-r)*(1-r)
    res2 = 1 + 6*r*r*(r-1)
    res = where(w>h, 0, where(w > 0.5 * h, 
        res1, res2))

    res *= 8/(pi*h*h*h)

    return res
