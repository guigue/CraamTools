class Phys(object): pass

class sb(Phys):
    v = 5.6705100000E-08
    d = 'Stefann Boltzmann'
    u = 'W m^-2 K-4'

class kb(Phys):
    v = 1.3806580000E-23
    d = 'Boltzmann'
    u = 'J K^-1'

class c(Phys):
    v = 2.9979245800E+08
    d = 'speed of light in vacuum'
    u = 'm s^-1'

class h(Phys):
    v = 6.6260755000E-34
    d = 'Planck'
    u = 'J s'

class G(Phys):
    v = 6.6725900000E-11
    d = 'Gravitational'
    u = 'm^2 Kg^-1 s^-2'

class e(Phys): pass

class mass(e):
    v = 9.1093897000E-31 
    d = 'electron rest mass'
    u = 'Kg'

    def keV(e):
        return 0.5109990615E+03

    class charge(e):
        v = 1.6021765650E-19
        d = 'electron charge'
        u = 'C'

class p(Phys): pass
class mass(p):
    v = 1.6726231100E-27
    d = 'proton mass'
    u = 'Kg'

class n(Phys): pass
class mass(n):
    v = 1.6749286100E-27
    d = 'neutron mass'
    u = 'Kg'

class Avo(Phys):
    v = 6.0221367000E+23
    d = 'Avogadro number'
    u = ''

class Conv(Phys): pass

class eV2K(Conv):
    v = 1.1604450000E+04
    d = 'eV to K convertion factor'
    u = 'K eV^-1'

class eV2J(Conv):
    v = 1.6021773300E-19
    d = 'eV to J'
    u = 'J eV^-1'

class Ang2keV(Conv):
    v = c.v * h.v / eV2J.v * 1.0E+07
    d = 'Angstrons to kilo eV'
    u = 'keV Angstrons^-1'
