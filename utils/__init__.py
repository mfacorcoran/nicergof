def nicer_plot_gti(gtifile, mkffile, parameter):
    """
    plots the gti intervals and MKF parameter for a given MKF parameter
    :parameter gtifile: name of the file containing the 'GTI' EXTENSION
    :parameter mkffile: name of the mkf file
    :parameter parameter: name of parameter in the mkf file to plot
    """
    from astropy.table import Table
    import matplotlib.pyplot as plt
    import numpy as np
    gti = Table.read(gtifile, hdu='GTI')
    gti['duration'] = gti['STOP'] - gti['START']
    mkf = Table.read(mkffile, hdu=1)
    t0 = gti['START'].min()
    fig = plt.figure(figsize=[10, 8])
    plt.xlim(-10, gti['STOP'].max() - t0)
    ymin = min(mkf[parameter])
    ymax = max(mkf[parameter])
    if ymin==ymax:
        ymin=0.0
        ymax=ymax*1.10
    plt.ylim(ymin, ymax)
    plt.ylabel(parameter.replace('_', '\_'), fontsize=16)
    plt.xlabel(f'Time - {t0:.3f}', fontsize=16)
    plt.step(mkf['TIME'] - t0, mkf[parameter], color='green')
    for b, e in zip(gti['START'], gti['STOP']):
        f = plt.fill(np.asarray([b, e, e, b, b]) - t0, [ymin, ymin, ymax, ymax, ymin], 'r', alpha=0.3)
    return

