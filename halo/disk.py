"""
Various potential formulations for disks (Miyamoto-Nagai, exponential)



"""
from numpy import power, square, pi

def mn_dens(M, a, r):
    """ 
    Surface density of the Miyamoto-Nagai disk when the scale height (b) is zero
    M - mass
    a - scale length
    r - (x^2 + y^2)^1/2

    returns Sigma (in units of M/r^2)
    """
    return M * a / ((a*a+r*r)**1.5 * 2*pi)

def mn_grav3(pos, M, a, b, G=43018.7):
    """ Gravity for Miyamoto-Nagai (disk) Potential """
    R2 = square(pos[:,:2]).sum(1)
    z2 = square(pos[:,2])
    acc = -G * M * reshape(power(R2+(a+(b*b+z2)**0.5)**2, -1.5), (len(pos),1))
    acc = acc*pos
    acc[:,2] *=1+a/(b*b+z2)**0.5

    return acc

def mn_grav(R, GM, a, b):
    """ Gravity for Miyamoto-Nagai (disk) Potential in plane"""
    return -R* GM * power(R*R+(a+b)**2, -1.5)

def mn_mean_dens(R, M, sum_ab):
    """ 
    Implied mean density from gravity in mn_grav, i.e. spherically averaged via

    4 pi r^2 rho_mean(r) = d/dr [ r^2 g / G]
    """

    return 3*M * square(sum_ab)/(4*pi*(R*R+square(sum_ab))**2.5)
