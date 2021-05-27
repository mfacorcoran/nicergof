__name__    = 'nicertimeconv'
__author__  = 'Mike Corcoran'
__version__ = '1.00'


def nicertimeconv(time, informat='met', outformat='mjd', outscale='utc',
                  LEAPINIT=2, TIMEZERO = 0, MJDREFI = 56658.0,
                  MJDREFF =  0.000777592592592593):
    """
    converts the input time to an output time in another time format

    converts MET to MJD via this recipe from CM

    On Oct 25, 2018, at 12:19 PM, Markwardt, Craig B (GSFC-6620) wrote:
        MJD(TT) = (MJDREFI+MJDREFF) + (TIMEZERO+TIME)/86400
        MJD(UTC) = (MJDREFI) + (TIMEZERO+TIME+LEAPINIT=2)/86400

    MJDREFI = 56658.0  corresponds to  Time("2014-01-01T00:00:00").mjd, the NICER epoch
    LEAPINIT = 2 is the number of leap seconds since the NICER epoch
    TIMEZERO = -1 seems to be the current value since mid 2018

    then uses astropy.time to convert to another specified output format

    (see NICER event lists section in https://heasarc.gsfc.nasa.gov/docs/nicer/mission_guide/)

    On Oct 25, 2018, at 12:19 PM, Markwardt, Craig B (GSFC-6620) <craig.b.markwardt@nasa.gov> wrote:
        MJD(TT) = (MJDREFI+MJDREFF) + (TIMEZERO+TIME)/86400
        MJD(UTC) = (MJDREFI) + (TIMEZERO+TIME+LEAPINIT=2)/86400


    :parameter time: input time (string or float)
    :parameter informat: format for input time
    :parameter outformat: output format for time
    :return:
    """
    from astropy.time import Time
    import numpy as np
    time = np.asarray(time)

    if informat == 'met':
        # convert time to mjd
        # This mjd is in TT
        mjd_tt = Time(MJDREFI + MJDREFF + (TIMEZERO + time)/86400., format='mjd', scale='tt')
        if outformat.lower().strip() != 'met':
            outtime = Time(mjd_tt, format='mjd', scale=outscale)
            try:
                setattr(outtime,'format',outformat)
            except Exception as errmsg:
                print('Could not convert {0} to {1} ({2})'.format(time,outformat, errmsg))
                print('Specified format must be one of')
                print(', '.join(Time.FORMATS))
    if outformat == 'met':
        time = Time(time, format=informat, scale='utc')
        mjd_tt = time.tt
        outtime = (mjd_tt.mjd - (MJDREFI+MJDREFF) )*86400. - TIMEZERO
    return outtime

if __name__ == '__main__':

    #isot = '2019-06-12T12:42:30'
    isot = ['2019-06-12T12:42:30', '2019-06-12T12:44:30']
    met = nicertimeconv(isot, informat='isot', outformat='met')
    for i, m in zip(isot, met):
        print('{0} UTC corresponds to NICER MET {1}'.format(i, m))

    met = [0, 2]
    isot = nicertimeconv(met, informat='met', outformat='isot')
    for i, m in zip(isot, met):
        print('{0} UTC corresponds to NICER MET {1}'.format(i, m))

    isot = '2014-01-01T00:00:00'
    met = nicertimeconv(isot, informat='isot', outformat='met')
    print('{0} UTC corresponds to NICER MET {1}'.format(isot, met))
    print(type(met))

    met = 3
    isot = nicertimeconv(met, informat='met', outformat='isot')
    print('{0} UTC corresponds to NICER MET {1}'.format(i, m))
    print(type(isot))

