#########################################################
#
# myGoes : some scripts to deal with GOES data in Python
#
#
# 
def get(isoStart,isoEnd):

    # Get: uses sunpy to get goes data.
    # isoStart, isoEnd: 'YYYY-MM-DD hh:mm:ss'
    # Returns file (object) and timeseries (object)
    
    from sunpy import timeseries as ts
    from sunpy.net import Fido
    from sunpy.net import attrs as a

    try:
        result = Fido.search(a.Time(isoStart, isoEnd), a.Instrument("XRS"))
        file=Fido.fetch(result)
        goes = ts.TimeSeries(file,concatenate=True)
    except:
        print("\n\n Data not found for the interval {0:s}--{1:s}".format(isoStart,isoEnd))
        file=[]
        goes = []

    return file, goes


