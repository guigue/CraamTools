#########################################################
#
# Goes : some scripts to deal with GOES data in Python
#        returns the filenames of raw data, and a dataframe with the data
#
#
# @Guiguesp - 2025-02-14BST19:41

__Version__ = '2025-02-14T19:41'

def get(isoStart,isoEnd):

    # Get: uses sunpy to get goes data.
    # isoStart, isoEnd: 'YYYY-MM-DD hh:mm:ss'
    # Returns file (object) and DataFrame (object)
    
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

    return file, goes.to_dataframe()


