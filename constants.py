UNITS = "SI"

UNIT_LENGTH = 1
UNIT_MASS   = 1
UNIT_TIME   = 1

DEFAULT_GRAVITATIONAL_CONSTANT  = 6.673e-11      # m3 kg-1 s-2          
DEFAULT_SPEED_OF_LIGHT          = 299792458      # m s-1
DEFAULT_SOLAR_MASS              = 1.98892e30     # kg
DEFAULT_PARSEC                  = 3.08568025e16  # m
DEFAULT_YEAR                    = 31556926       # s

DEFAULT_h = 0.73

G = GRAVITATIONAL_CONSTANT = DEFAULT_GRAVITATIONAL_CONSTANT
c = SPEED_OF_LIGHT = DEFAULT_SPEED_OF_LIGHT
h = DEFAULT_h

def set_h(new_h):
    global h
    h = new_h
    set_units(UNITS)

def set_units(units):
    global UNITS
    global c, SPEED_OF_LIGHT, G, GRAVITATIONAL_CONSTANT
    
    if units=="SI":
        UNIT_LENGTH = 1
        UNIT_MASS = 1
        UNIT_TIME = 1

    elif units=="GALACTIC":
        UNIT_LENGTH = (1e6 * DEFAULT_PARSEC / h)                # 1.0 Mpc h^-1
        UNIT_MASS = (1e10 * DEFAULT_SOLAR_MASS / h)             # 10^10 M_solar h^-1
        UNIT_TIME = (1e3 * DEFAULT_PARSEC / h)                  # 977.8 Gyr h^-1
    
    elif units=="CGI":
        UNIT_LENGTH = 0.01
        UNIT_MASS = 0.001
        UNIT_TIME = 1

    UNITS = units
    G = GRAVITATIONAL_CONSTANT = DEFAULT_GRAVITATIONAL_CONSTANT * UNIT_MASS * UNIT_TIME**2 / UNIT_LENGTH**3
    c = SPEED_OF_LIGHT = DEFAULT_SPEED_OF_LIGHT * UNIT_TIME / UNIT_LENGTH;
        
set_units("SI")