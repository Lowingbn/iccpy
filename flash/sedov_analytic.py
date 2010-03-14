"""
Peter Creasey
p.e.creasey.00@googlemail.com

solution to the Sedov problem

based on the C code by Aamer Haque
"""
from scipy.special import gamma as Gamma
from numpy import power, arange, empty, float64, log, exp, pi, diff, inner
import pylab as pl

def calc_a(g=1.4,nu=3):
    """ calc some constants for the sedov solution """
    a = [0]*8
   
    a[2] = (1.0-g) / (2.0*(g-1) + nu)
    a[3] = nu / (2.0*(g-1.0) + nu)
    a[5] = 2.0 / (g-2.0)
    a[6] = g / (2.0*(g-1.0) + nu)
   
    a[1] = ( ((nu+2.0)*g)/(2.0+nu*(g-1.0)) ) * ( (2.0*nu*(2.0-g))/(g*(nu+2.0)**2) - a[2])
    a[4] = (a[1]*(nu+2.0)) / (2.0-g)
    a[7] = (2.0 + nu*(g-1))*a[1]/(nu*(2.0-g))
    return a

def calc_beta(v, g=1.4, nu=3):
    """ beta values for the sedov solution (exponents of the polynomials) """
    beta = empty((4,v.size), float64)
    beta[0] = (nu+2.0)*(g+1.0)*v/4.0
    beta[1] = ((g+1.0)/(g-1.0)) * ( (nu+2.0)*g*v/2.0 - 1.0 )
    beta[2] = ((nu+2.0)*(g+1.0)) / ( (nu+2.0)*(g+1.0) -2.0*(2.0 + nu*(g-1.0)) )
    beta[2] *= 1.0 - (2.0 + nu*(g-1.0))*v/2.0
    beta[3] = ((g+1.0)/(g-1.0)) * (1.0 - (nu+2.0)*v/2.0)
    return beta


def sedov(t, E0, rho0, g, n=10000, nu=3):
    """ 
    solve the sedov problem
    t - the time
    E0 - the initial energy
    rho0 - the initial density
    n - number of points (10000)
    nu - the dimension
    g - the polytropic gas gamma
    """
    # the similarity variable
    v_min = 2.0 / ((nu + 2.0) * g)
    v_max = 4.0 / ((nu + 2.0) * (g + 1.0))

    v = v_min + arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = log(beta)
    
    r = exp(-2.0 / (2.0 + nu) * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * exp(a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2])
    p = exp((2 * nu / (2.0 + nu)) * lbeta[0] + (a[5] + 1) * lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2])
    u = beta[0] * r * 4.0 / ((g + 1) * (nu + 2))
    p *= 8.0 / ((g + 1) * (nu + 2) * (nu + 2))

    # volume of an n-sphere
    dv = pi ** (nu / 2.0) * power(r, nu) / Gamma(nu / 2.0 + 1)
    dv[1:] = diff(dv)
    dv[0] = 0.0

    # (dimensionless) energy of the model solution
    q = inner(rho * u * u * 0.5 + p / (g - 1), dv)

    # the factor to convert to this particular problem
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))

    # solve for the shock speed
    shock_speed = fac * (2.0 / (nu + 2))
    rho_s = ((g + 1.0) / (g - 1.0)) * rho0                                                                            
    r_s = shock_speed * t * (nu + 2) / 2.0
    p_s = (2.0 * rho0 * shock_speed ** 2) / (g + 1.0)
    u_s = (2.0 * shock_speed) / (g + 1.0)
    
    r *= fac * t
    u *= fac
    p *= fac * fac * rho0
    rho *= rho0
    return (r, p, rho, u), (r_s, p_s, rho_s, u_s, shock_speed)


def test():
    """ draw a 3d sedov solution """
    val, sval =  sedov(t=0.05, E0=5.0, rho0=5.0, g=5.0/3.0)
    r,p,rho,u = val
    r_s, p_s, rho_s, u_s, shock_speed = sval
    print 'rho shock', rho_s
    print 'p shock', p_s
    print 'u shock', u_s
    print 'r shock', r_s
  
    dv = (4/3.0)*pi*r*r*r
    dv[1:] = diff(dv)

    # thermal and kinetic energy
    te = (p*dv/(5.0/3.0-1))
    ke = (rho*u*u*0.5*dv)
    energy = te.sum() + ke.sum()
    mass = inner(rho,dv)
    print 'density', mass / (4/3.0 * pi * r_s**3)
    print 'energy', energy

    print 'shock speed', shock_speed
    pl.plot(r/r_s,rho/rho_s, label=r'$\rho/\rho_s$')
    pl.plot(r/r_s,p/p_s, label=r'$p/p_s$')
    pl.plot(r/r_s,u/u_s, label=r'$u/u_s$')
    pl.legend(loc='upper left')
    pl.show()
        
if __name__=='__main__':
    test()