from astropy.table import Table

def lc2gti(lcname):
    """
    returns a Table object from a HEASARC-formatted RATE file
    with columns START, STOP and DURATION added
    :param lcname: name of HEASARC-formatted FITS RATE file
    :return: table object of the RATE extension with START, STOP and DURATION added
    """
    lctab = Table.read('ni2012080116.lc', hdu='RATE')
    lctab['START'] = lctab['TIME'] + lctab.meta['TIMEZERO'] - lctab.meta['TIMEDEL'] * lctab.meta['TIMEPIXR']
    lctab['STOP'] = lctab['TIME'] + lctab.meta['TIMEZERO'] + lctab.meta['TIMEDEL'] * lctab.meta['TIMEPIXR'] * lctab[
        'FRACEXP']
    lctab['DURATION'] = lctab['STOP'] - lctab['START']
    return lctab