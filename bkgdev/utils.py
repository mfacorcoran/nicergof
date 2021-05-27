"""
so here's what we do:
1) get the good times from the spectrum
2) read in the bkg events file as a Table & make mkftab, table of the mkf3 file for the observation
3) filter the mkf3 table for the first good time (START <= TIME <= STOP) - use mkf_tfilt
    3.1) from the time-filtered mkf file, get the KP distribution
    3.2) for the first KP bin (KPlo, KPhi) with data
        3.2.1) select the events from the bkg evt table with KPlo<=KP<=KPhi as a table (tabsel)
        3.2.2) get the COR_SAX distribution for the first KP bin
        3.2.3) for the first COR_SAX bin
            3.2.2.1) select the events from tabsel with CORlo<=COR_SAX<=COR_HI as a table
            3.2.2.2)
"""

import numpy as np

def mkf_tfilt(mkffile, tstart, tstop, extname = 'PREFILTER' ):
    """
    filters an mkf file so that tstart<=MKF['time'] <= tstop
    Useful for selecting MKF rows within a GTI, for example
    """
    from astropy.io import fits
    from astropy.table import Table
    import numpy as np
    mkf = fits.open(mkffile)
    mkfdata = mkf[extname].data
    mkftime = mkfdata.TIME
    it = np.where((mkftime >= tstart) & (mkftime <= tstop))[0]
    if it.size == 0:
        print("Could not find MKF times between {0:.3f} and {1:.3f} in {2}".format(tstart, tstop, mkffile))
        return None
    mkf[1].data = mkf[1].data[it]
    mkftab = Table(mkf[1].data)
    return mkftab

def tabselect(tab, colname, lo, hi):
    """
    filters a table on column name colname, selecting those events where lo<=(colname value)<=hi
    :param tab:  table (astrop.table.Table)
    :param colname: column name of table to select on
    :param lo: lower limit
    :param hi: upper limits
    """
    tabsel = tab[(tab[colname]>=lo) & (tab[colname]<=hi)]
    return tabsel

def get_tdurs(time, deltat = 1.0, verbose=False):
    """
    calculates durations of discrete time intervals for an array of times
    :param time: array of times taken at non-continuous intervals
    :param deltat: time bin size which defines a "continuous  interval". By default continuous intervals are < 1 second in length
    :param verbose: if True print diagnostic messages
    """
    #
    #
    #
    a = list(time)
    a.sort()
    a=np.asarray(a)
    ad = a[1:]-a[:-1]
    try:
        ints = np.where(ad>deltat)[0]
    except:
        # if where didn't find any values, i.e. all time values occur within deltat
        return max(time)-min(time)
    iints = np.asarray(ints)+1
    # initial interval
    tdur = a[iints[0]-1]-a[0]
    print(f"initial duration: {tdur:.2f} seconds")
    for i in enumerate(iints[:-1]):
        i=i[0]
        #print(iints[i], iints[i+1])
        dt = a[iints[i+1]-1]- a[iints[i]]
        tdur += dt
        if verbose:
            #print(dt, a[iints[i]:iints[i+1]])
            print(f" total duration now {tdur:.2f}")
    # final interval
    td = a[-1]-a[iints[-1]]
    tdur += td
    print(f"Total duration now {tdur:.2f}")
    return tdur
