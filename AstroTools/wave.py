import CraamTools.AstroTools.constants as cte
import pdb

class wave:

    def __init__(self,val=30,unit="THz"):
        
        if (unit == "angstrom"):
            _fac = 1.0E-10
            _type = "wl"
        elif (unit == "nm"):
            _type = "wl"
            _fac = 1.0E-09
        elif (unit == "um"):
            _type = "wl"
            _fac = 1.0E-06
        elif (unit == "mm"):
            _type = "wl"
            _fac = 1.0E-03
        elif (unit == "cm"):
            _type = "wl"
            _fac = 1.0E-02
        elif (unit == "GHz" or unit == "ghz"):
            _type = "fq"
            _fac = 1.0E+09
        elif (unit == "THz" or unit == "thz"):
            _type = "fq"
            _fac = 1.0E+12
        elif (unit == "MHz" or unit == "mhz"):
            _type = "fq"
            _fac = 1.0E+06
        else:
            _type = "wl"
            _fac = 1.0

        self.unit = "SI"
        if (_type=="wl"):
            self.wl = val * _fac
            self.fq = cte.Phys.c / (val * _fac)

        elif (_type=="fq"):
            self.fq = val * _fac
            self.wl = cte.Phys.c / (val* _fac)
            
        return
