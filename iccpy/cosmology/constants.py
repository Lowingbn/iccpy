import numpy as np
import iccpy.constants
import scipy.integrate

OMEGA_L0 = OmegaLambda0            = 0.742
OMEGA_M0 = OmegaMatter0            = 0.238
OMEGA_B0 = OmegaBaryon0            = 0.0418
OMEGA_K0                           = 0
OMEGA_R0                           = 0
h                                  = 0.73

_times = {}

def _hubble_func(a):
    return np.sqrt(OMEGA_K0/(a*a) + OMEGA_M0/(a**3) + OMEGA_R0/(a**4) + OMEGA_L0);

def hubble_0(h=0.73):
    return h * 3.2407764868e-18 * iccpy.constants.UNIT_TIME
    
def hubble_param(a, h=0.73):
    return hubble_0(h) * _hubble_func(a)

def time(a):
    if a not in _times: 
        _times[a] = scipy.integrate.quad(lambda a: 1.0/(a*_hubble_func(a)), 0, a, ())[0]/hubble_0()
    
    return _times[a]

def omega_m(a):
    h = _hubble_func(a)
    return OMEGA_M0/(a**3 * h**2)
    
def omega_l(a):
    h = _hubble_func(a)
    return OMEGA_L0/(h**2)

def omega_r(a):
    h = _hubble_func(a)
    return OMEGA_R0/(a**4 * h**2)

def omega_k(a):
    h = _hubble_func(a)
    return OMEGA_K0/(a**2 * h**2)

def omega_crit(a):
    H = hubble_param(a)
    return 3 * H**2 / (8 * np.pi * iccpy.constants.G) 

def omega_crit0():
    return 3 * hubble_0()**2 / (8 * np.pi * iccpy.constants.G) 
