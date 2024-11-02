###########################
#
# Constants
#
# Taken from "Allen's Astrophysical Quantities",
# 4th Edition. Edited by Arthur N. Cox.
# ISBN: 0-387-98746-0, published 1999
#
# All constants are in SI units except especifically
# indicated.
#
# @ guiguesp - Sampa , 2018-03-04 
#
##########################

class Phys:
#
# Constants taken from NIST site. 
# 2019-05-20
#
    Stefann_Boltzmann = 5.6705100000E-08
    Boltzmann         = 1.3806490000E-23  
    c                 = 2.9979245800E+08  # speed of light
    Planck            = 6.6260701500E-34
    G                 = 6.6725900000E-11  # Newton gravitational constant
    e_mass            = 9.1093897000E-31  # electron mass
    e_rmass           = 0.5109990615E+03  # electron rest mass (keV)
    e_charge          = 1.6021766340E-19  # elementary charge
    p_mass            = 1.6726231100E-27  # proton mass
    n_mass            = 1.6749286100E-27  # neutron mass
    Avogadro          = 6.0221407600E+23  # Avogadro's Number
    eV2K              = 1.1604450000E+04  # eV --> Kelvin
    ev2J              = 1.6021773300E-19  # eV --> Joule
    A2keV             = c*Planck/ev2J/1.0E-04
    
class Astro:
    
    ly           = 9.4607304720E+15 # light year
    pc           = 3.0856776000E+16 # parsec
    au           = 1.4959787066E+11 # Astronomical Unit
    s_mass       = 1.9891000000E+30 # Sun mass
    s_radius     = 6.9550800000E+08 # Sun Radius
    s_luminosity = 3.8450000000E+26 # Sun luminosity


