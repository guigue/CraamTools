# cd $HOME/solar/2017
import oRBD
d = oRBD.RBD()
d.readRBDinDictionary('rs1170819.1400')

from CraamTools.AstroTools import sstMap
import importlib

importlib.reload(sstMap)
m0=sstMap.sstMap(d)
m0.Scans.Extract_Scans(d,0)
m0.Scans.Correct_Scans(eshift=5)
m0.Get_Observed_Sun_Center()
m0.Map.Create_Image(m0.Scans,m0.Observed_Sun_Center)
m0.Map.Rotate_North_Image(m0.PB0,m0.parallactic_angle)
m0.Calibrate_Map_in_T()
