# astronomical units in cgs
from math import pi, e

c = 2.9979245800000e10 # speed of light, cm/s
pc = 3.0856775807000e18 # 1 parsec, cm
yr = 3.1556925200000e07 # year in s

m_p = 1.6726231100000e-24 # proton mass, g
G = 6.6725985000000e-08 # big G, cm^3 s^-2 g^-1
k_b = 1.3806581200000e-16 # boltzmann constant cm^2 s^-2 g K^-1 
 
sigma = 5.6705119000000e-5 # stefan boltzmann constant,  s^-3 g K^-4 
msun = 1.9889225000000e33 # solar mass, g  

h = 6.6260755400000E-27 # planck constant cm^2 s^-1 g 

kpc = 1e3 * pc
mpc = 1e3 * kpc

myr = 1e6 * yr


##-----------------------------------------------------------
# Cosmological Parameters


HubbleConstant         = 2.3657686e-18   # 0.73
OmegaMatter            = 0.238
OmegaBaryon            = 0.0418

rho_crit_z0 = 3 * HubbleConstant**2 / (8 * pi * G)



