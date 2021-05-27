def combine_mkf(mkflist, outfile, clobber=True):
    """
    For a list of nicer mkffiles, combine the data, time sort it, and write to outfile
    :param mkflist: input list of mkf files
    :param outfile: output combined mkffiles
    :return: astropy Table   of sorted, combined mkf data
    """
    from astropy.table import Table, vstack
    tablist = []
    for m in mkflist:
        tablist.append(Table.read(m, hdu='PREFILTER'))
    combinedtab = vstack(tablist)
    combinedtab.sort(['TIME'])
    combinedtab.write(outfile, format='fits', overwrite=clobber)
    return combinedtab
