def check_gain(obsid, rootdir='.'):
    """
    checks if an obsid uses the most recent gain file
    """
    from astropy.io import fits
    import os
    import glob
    from subprocess import run
    ufafile = glob.glob(os.path.join(rootdir,str(obsid),'xti', 'event_cl','ni{0}_0mpu7_ufa.evt*'.format(obsid)))[0]
    try:
        gainused = fits.open(ufafile)[1].header['GCALFILE']
    except Exception as e:
        print('Could not open {0} ({1}); returning'.format(ufafile, e))
        return -1
    print('Gain file used: {0}'.format(gainused))
    curgfile = run(['quzcif', 'nicer', 'xti', '-', '-', 'MPU_GAIN', 'now', 'now', '-'], capture_output=True, text=True).stdout
    curgfile = curgfile.split()[0].split('/')[-1]
    print('Latest Gain File: {0}'.format(curgfile))
    if gainused == curgfile:
        print('Events processed with latest gain file')
    else:
        print('Events NOT processed with latest gain file')
        print('See https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/gain-cal/ for details')
    return 0