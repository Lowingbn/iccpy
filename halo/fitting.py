import numpy as np

def density_broken_power_law(r, rho0, r_br, n_inner, n_outer):
    rho1 = rho0 * r_br**n_inner/r_br**n_outer
    rho = np.where(r<r_br, rho0*r**n_inner, rho1*r**n_outer)
    return rho

def density_power_law(r, rho0, n):
    rho0 = np.max([1e-10, rho0])
    return rho0 * r**n

def leastsq_power_law_residuals(p, den, r, dV):
    rho = density_power_law(r, *p)

    error = np.sqrt(rho/dV)
    return (rho-den)/error

def leastsq_broken_power_law_residuals(p, den, r, dV):
    rho = density_broken_power_law(r, *p)

    error = np.sqrt(rho/dV)
    return (rho-den)/error

def fit_power_law(density, r, dV, p0=None):
    from scipy.optimize import leastsq

    if p0 is None:
        p0 = (1.0, 1.0)
    p, plsq = leastsq(leastsq_power_law_residuals, p0, args=(density, r, dV), maxfev=10000)

    chi2 = np.sum(leastsq_power_law_residuals(p, density, r, dV)**2)
    return p[0], p[1], chi2

def fit_broken_power_law(density, r, dV, p0=None):
    from scipy.optimize import leastsq

    if p0 is None:
        p0 = (1.0, 1.0, 1.0, 1.0)
    p, plsq = leastsq(leastsq_broken_power_law_residuals, p0, args=(density, r, dV), maxfev=10000)

    chi2 = np.sum(leastsq_broken_power_law_residuals(p, density, r, dV)**2)
    return p[0], p[1], p[2], p[3], chi2

if __name__=="__main__":
    r = np.logspace(-2, 1, 100)
    x = np.where(r<1, r**2, r**3)

    print fit_power_law(x, r, 1)
    print fit_broken_power_law(x, r, 1)
