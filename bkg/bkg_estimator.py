"""
This package contains functions which can be used to create an estimated NICER background spectrum based on the "environmental" model
developed by the NICER Guest Observer Facility. The environmental model uses a combination of the cut-off rigidity (COR_SAX)
and the Planetary K index (KP) which gives an estimate of the space weather environment.  This model also uses the SUN_ANGLE parameter
which helps describe the low-energy background produced by optical loading.  COR_SAX and SUN_ANGLE are contained in the
"make filter" file (either the standard auxil/ni*.mkf file distributed with processed data or -recommended- the augmented MKF file produced by
the "niprefilter2" tool distributed with the NICERDAS HEASoft package).  The KP values are not currently included in either the .mkf or the
enhanced (.mkf2) makefilter files, and must be added using the add_kp function defined here.

PRELIMINARIES:

You'll need access to these files:
    a) The background events file 30nov18targskc_enhanced.evt (current version: https://heasarc.gsfc.nasa.gov/FTP/caldb/data/nicer/xti/pcf/30nov18targskc_enhanced.evt)
    b) the KP.fits file (the current version, updated daily, is at https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/kp.fits)
You can access these files virtually by specifying the URLs given above (the default for the functions defined below) or
you can download them to a local directory for faster access.

EXAMPLE of creating an estimated NICER background spectrum:

  Assumptions:
  - the NICER data is stored in subdirectories of a root directory called "testdata" of the current working directory.
  - the NICER observation obsid is 1200040103
  - a extracted spectrum (which contains background) is called src.pha and is located in the current working directory.
  - a NICER make filter file created using the niprefilter2 NICERDAS tool exists.

    1) Use the niprefilter2 tool from NICERDAS to create a .mkf2 file from the .mkf file in the <obsid>/auxil directory (where <obsid> is the nicer  10-digit observation ID, and the standard NICER data directory structure is assumed)

    2) after making sure that the nicergof/bkg directories are in your $PYTHONPATH, import this package as "be":

        >>> from nicergof.bkg import bkg_estimator as be

    3) update the .mkf2 file to include the KP values, using the "add_kp" function (assume this is for an observation with obsid=1200040103):

        >>> mkf3 = be.add_kp("testdata/1200040103/auxil/ni1200040103.mkf2")

    This will create a ".mkf3" file in the 11200040103/auxil directory

    4) use the mk_bkg_spec_evt function to create the background spectrum:

        >>> bkg_chan, bkgspectot, btotexpo = mk_bkg_spec_evt('src.pha', mkf3file=mkf3)

    This will create a HEASARC-compliant background PHA file (with .pha replaced by _bkg.pha, i.e. the background file for "./nicer.pha" is
    "./nicer_bkg.pha") and also return the background channels, spectrum and exposure as the python variables bkg_chan, bkgspectot, btotexpo.

CAVEATS:
    * This is PRE-RELEASE software.
    * This software estimates the instrumental background.  X-ray cosmic background (sky background) appropriate to your source is NOT included in the estimated background spectrum produced by mk_bkg_spec_evt().  However, cosmic X-ray background in the NICER blank fields is included in the estimated background.
    * Note that the background events file excludes FPMs #14 & 34, so it uses data from 50 out of the 52 active FPMs.  The estimated background
    * There may be combinations of (KP, COR_SAX, SUN_ANGLE) which are not contained in the background events file; in this case these times are ignored in the output background spectrum.
    * there may be other parameters that are important in determining background, or other parameters which give a better estimate of the background.  This is still under investigation.
    * there are undoubtedly other issues.

"""

__author__ = "M. F. Corcoran (NASA/GSFC & CUA)"
__version__ = "0.2"
__status__ = "Pre-release"

import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.time import Time
import os
import pandas as pd
import glob


def add_kp(mkffile, kpfile ='https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/kp.fits',
           extname = 'PREFILTER', outname = '', clobber=False, verbose=True):
    """
    updates an mkf file for the kp values from the kp.fits file
    writes out the file with an extension outexten (".mkf3" by default)
    :param mkffile: name of .mkf2 file
    :param kpfile: name of FITS file with KP vs TIME
    :param clobber: if true overwrite existing file

    :return: writes a ".mkf3" file and returns a status (0 if success)
    """
    status = 0
    kphdu = fits.open(kpfile)
    kp = kphdu[1].data.KP
    kptime = kphdu[1].data.TIME
    mkf = fits.open(mkffile)
    mkfdata = mkf[extname].data
    mkfhead = mkf[extname].header
    mkftime = mkfdata.TIME
    mkfcols = mkfdata.columns
    kpinterp = np.interp(mkftime, kptime, kp)
    if 'KP' in mkf[extname].data.columns.names:
        # if KP column already exists in mkf file just update the column
        mkupdatedhdu = mkf
        mkupdatedhdu[extname].data['KP'] = kpinterp
    else:
        # append column to mkf hdu
        # see http://docs.astropy.org/en/stable/io/fits/usage/table.html "Merging Tables"
        # define a new KP column
        kpcol = fits.ColDefs([fits.Column(name='KP', format='D', array=kpinterp)])
        mkupdatedhdu = fits.BinTableHDU.from_columns(mkfcols + kpcol, header=mkfhead)
    mkfdir = os.path.split(mkffile)[0]
    if len(outname) == 0:
        mkfoutfile = os.path.split(mkffile)[-1].split('.')[0]+'.mkf3'
        outname = os.path.join(mkfdir, mkfoutfile)
    if verbose:
        print("Writing {out}".format(out=outname))
    try:
        mkupdatedhdu.writeto(outname, overwrite=clobber)
    except Exception as e:
        print(e)
        status = -1
        pass
    return status


def mk_bkg_spec_evt(srcpha, mkf3file, bevt="https://heasarc.gsfc.nasa.gov/FTP/caldb/data/nicer/xti/pcf/30nov18targskc_enhanced.evt",
                    numfpms = 50, gtimin=0, gti=None, clobber=True, verbose=True):
    """
    This routine creates a NICER instrumentental background spectrum from the NICER background event file (provided by K. Gendreau) for a specified extracted spectrum (containing source + background counts).  It does this by looking at the range of KP, SUN_ANGLE, COR_SAX for the observation, then creating histogram bins of these values, extracting events from the background event file which match the observed KP, SUN_ANGLE and COR_SAX bin range, and correcting the background exposure for dead time and amount of time the background observation was in that particular (KP, SUN_ANGLE, COR_SAX) bin, and finally correcting for the number of Focal Plane Modules (FPM) in use

    :param srcpha: name of the pha file of the source spectrum (with full path, with an extension of ".pha")
    :param mkf3file: the name of the MKF file for the observation
    :param bevt: the name of the background event file
    :param numfpms: number of FPMS from which the NICER source spectrum was extracted.  Typically 50 out of 52 active FPMs are used in data extraction (FPM 14 and 34 are usually excluded since they tend to be noisier).
    :param gtimin: minimum acceptable time in a GTI in seconds.  GTIs less than this value will not be included in the calculated bkg
    :param gti: good time interval table from the observation as a table object or dataframe with at minimum columns = ['START', 'STOP'], used to select background from a particular time interval
    :return: bkg_chan, bkgspec_tot, btotexpo_kcsa; bkg_chan is the NICER channel array corresponding to the total background
     spectrum (bkgspec_tot), and btotexpo_kcsa (total exposure in the bkg events after filtering on KP, COR & SA)
     Writes out the background PHA file to a file in the same directory as the pha file with the extension of the pha file replaced by "_bkg.pha"
    """
    # read the "enhanced" background events table (includes kp, cor, and sun angle)
    betab = Table.read(bevt, hdu=1)
    # and get total bkg exposure
    bkgexpotot = fits.open(bevt)[1].header['EXPOSURE']

    if gti is None:
        gti = fits.open(srcpha)['GTI'].data
        gti = pd.DataFrame(gti)
        gti['Duration'] = gti.STOP - gti.START

    # GET THE MKF3 INFO FOR THIS OBSERVATIONS
    # mkf3 file from the observation
    mkf3 = fits.open(mkf3file)
    # these are the times, kp, SA, cor and deadtime values for the GTI from the observation
    try:
        mkfkp = mkf3['PREFILTER'].data.KP
    except:
        print("Could not get KP information; returning")
        return 0, 0, 0
    mkftime = mkf3['PREFILTER'].data.TIME
    mkfSA = mkf3['PREFILTER'].data.SUN_ANGLE
    mkfcor = mkf3['PREFILTER'].data.COR_SAX
    mkfDT = mkf3['PREFILTER'].data.MPU_DEADTIME.mean(axis=1)

    # for a given GTI, find the start and stop times for the gti, then get the values of
    # KP, COR_SAX, SUN_Angle, deadtime within the GTI

    for gtinum in range(len(gti)):
        skipit = False
        gdur = gti.iloc[gtinum].Duration
        print("\nFor GTI #{gn}; Duration = {d}\n".format(gn=gtinum, d=gdur))
        if gdur > gtimin:

            t0s = gti.iloc[gtinum].START
            t0e = gti.iloc[gtinum].STOP

            # it are the indices which select the specified good time interval
            it = np.where((mkftime >= t0s) & (mkftime <= t0e))[0]
            if it.size == 0:
                print("Could not find MKF times between {0:.3f} and {1:.3f} in {2}".format(t0s, t0e, mkf3file))
                skipit = True
            if not skipit:
                # now find time length of time that kp > kplo and kp < kphi
                gtidur = gti.loc[gtinum].Duration

                # select parameter values (kp, cor, sa, dt) corresponding to the selected GTI
                mkfkpgti = mkf3['PREFILTER'].data.KP[it]
                mkftimegti = mkf3['PREFILTER'].data.TIME[it]
                mkfSAgti = mkf3['PREFILTER'].data.SUN_ANGLE[it]
                mkfcorgti = mkf3['PREFILTER'].data.COR_SAX[it]
                mkfDTgti = mkf3['PREFILTER'].data.MPU_DEADTIME.mean(axis=1)[it]

                # NOW SELECT BY KP values for the given gtinum

                a = np.histogram(mkfkpgti, bins=np.arange(0, 15))
                kpnum = a[0]
                kpwt = kpnum / kpnum.sum()  # kpwt is the fraction of the gti duration in that KP bin
                kpbin = a[1]

                # for each non-zero kp bin get the fraction of the gti duration
                # that the kp value was within that KP bin

                # ik selects only those kpwt bins which are non-zero
                ikk = np.where(kpwt > 0.0)[0]

                for indk in ikk:
                    #  kpbin[i] gives the lower edge of the kp bin for the i'th (non-empty) bin
                    #  kplo, kphi give the lower, upper bin edge
                    #  for example ik[0] is the index in the kpbin of the first non-zero bin
                    kplo = kpbin[indk]
                    kphi = kpbin[indk + 1]

                    # since the mkf data are sampled once per second, the total time in this KP bin is just the
                    # length of the isel array * 1 second.  Then we have to multiply that time by the
                    # fraction of time the kp value was in that bin, i.e. by kpwt[i]
                    isel = np.where((mkfkpgti >= kplo) & (mkfkpgti <= kphi))[0]

                    # select parameter values (kp, cor, sa, dt) corresponding to the kp bin range
                    mkfkpsel = mkf3['PREFILTER'].data.KP[isel]
                    mkftimesel = mkf3['PREFILTER'].data.TIME[isel]
                    mkfsasel = mkf3['PREFILTER'].data.SUN_ANGLE[isel]
                    mkfcorsel = mkf3['PREFILTER'].data.COR_SAX[isel]
                    mkfdtsel = mkf3['PREFILTER'].data.MPU_DEADTIME.mean(axis=1)[isel]

                    # for now let's assume the times are continuous;
                    #  TODO: check for discontinuities
                    # to get the total exposure in this kp bin, get the total number of seconds = isel*0.0+1 (first term),
                    #   correct for Deadtime (second term),
                    #   and multiply by the fraction of total time in the selected kp bin (third term)
                    expo = np.sum((isel * 0.0 + 1.0) * (1.0 - mkfDT[it[isel]]) * kpwt[indk])
                    # if the exposure is more than gtidur, truncate it to gtidur
                    # note: this should only be needed when the exposure is just a bit larger than gtidur; if it's much larger then there
                    # are problems!
                    if expo > gtidur:
                        expo = gtidur
                    if verbose:
                        print("KP: expo = {ex}; dur of GTI = {gtid}".format(ex=expo, gtid=gtidur))

                    # so we've extracted the times when the KPs are within the first bin range;
                    # we NEXT need to extract BKG events from the bkg events that have KPs in the specified range,
                    # then accumulate the events and the exposure time -
                    # (assume a single "GTI" of duration equal to the exposure time)

                    betabselk = betab[(betab['KP'] > kplo) & (betab['KP'] <= kphi)]  # number of bkg events
                    bhist = np.histogram(betab['KP'], bins=kpbin)[0]
                    bexpowt = bhist[indk] / float(bhist.sum())  # fraction of bkg exposure in this kp bin
                    bdtave = betabselk['MPU_DEADTIME'].mean()  # use the average deadtime for bkg exposure calculation
                    bexpototk = bkgexpotot * (1 - bdtave) * bexpowt

                    # then we want to look at the distribution of cor
                    # for all the cor values in the mkf3 for all values with kp within the selected kp bin,
                    #    find the distribution of cor values
                    #      then for each non-empty cor bin
                    #        * calculate fraction of expo of time where cor is in that bin
                    #        * extract bkg events from the kp-selected events which have cor values within that bin

                    # print(bexpototk, bkgexpotot)
                    if verbose:
                        print("KP: % of bkg events between {kpl} < kp < {kph} = {bper:.2f}%".format(kpl=kplo, kph=kphi,
                                                                                                bper=bexpototk / bkgexpotot * 100))

                    #
                    # NOW FOR COR; bexpotot is the total background exposure in the kp bin - now we need the fraction of that exposure
                    # which has COR in the selected COR bin
                    #

                    # set the initial exposure to the KP bin exposure
                    cexpo = expo

                    corbins = np.arange(1, 20)
                    a = np.histogram(mkfcorsel, bins=corbins)
                    cornum = a[0]
                    corwt = cornum / float(cornum.sum())  # corwt is the fraction of the gti duration in that cor bin
                    corbin = a[1]

                    # for each non-zero kp bin get the fraction of the gti duration
                    # that the kp value was within that KP bin

                    # ik selects only those kpwt bins which are non-zero
                    ikc = np.where(corwt > 0.0)[0]

                    for indc in ikc:
                        # bin ranges for this bin
                        corlo = corbin[indc]
                        corhi = corbin[indc + 1]

                        # since the mkf data are sampled once per second, the total time is just the
                        # length of the isel array * 1 second.  Then we have to multiply that time by the
                        # fraction of time the kp value was in that bin, i.e. by kpwt[i]
                        isel = np.where((mkfcorsel >= corlo) & (mkfcorsel <= corhi))[0]

                        # for now let's assume the times are continuous; TODO:check for discontinuities
                        # include Deadtime correction (second term) and fraction of total time in kp,cor bin (third term)
                        cexpo = ((isel * 0.0 + 1) * (1.0 - mkfdtsel[isel]) * corwt[indc]).sum()
                        if cexpo > gtidur:
                            cexpo = gtidur
                        if verbose:
                            print("COR: expo = {ex}; dur of GTI = {gtid}".format(ex=cexpo, gtid=gtidur))

                        # so we've extracted the times when the KPs are within the first bin range;
                        # and the cor's are in their selected bin
                        # we need to extract BKG events from the bkg events that have KPs and CORs in the specified bins,
                        # then accumulate the bkg events and the bkg exposure time -
                        # (assume a single bkg "GTI" of duration equal to the exposure time)

                        # number of bkg events in the selected kp, cor bin
                        betabselc = betabselk[(betabselk['COR_SAX'] > corlo) & (betabselk['COR_SAX'] <= corhi)]
                        bhist = np.histogram(betabselk['COR_SAX'], bins=corbin)[0]
                        bexpowt = bhist[indc] / float(bhist.sum())  # fraction of bkg exposure in this kp, cor bin
                        bdtave = betabselc['MPU_DEADTIME'].mean()  # use the average deadtime for bkg exposure calculation
                        bexpototc = bexpototk * (1.0 - bdtave) * bexpowt

                        # then we want to look at the distribution of cor
                        # for all the cor values in the mkf3 for all values with kp within the selected kp bin,
                        #    find the distribution of cor values
                        #      then for each non-empty cor bin
                        #        * calculate fraction of expo of time where cor is in that bin
                        #        * extract bkg events from the kp-selected events which have cor values within that bin

                        # print(bexpototc, bexpototk, bkgexpotot)
                        if verbose:
                            print("COR: % of bkg events between {bl} < COR < {bh} = {bper:.2f}%".format(bl=corlo, bh=corhi,
                                                                                                    bper=bexpototc / bkgexpotot * 100))
                            print("Done COR index {0}; {1} < COR < {2}".format(indc, corlo, corhi))

                        #
                        # NOW FOR SUN ANGLE
                        #

                        sexpo = cexpo

                        sabin = np.arange(40, 190, 10)
                        a = np.histogram(mkfsasel, bins=sabin)
                        # q = hist(mkfsasel, bins = sabin)
                        sanum = a[0]
                        sawt = sanum / float(sanum.sum())  # the fraction of the gti duration in that SA bin
                        sabin = a[1]

                        # for each non-zero kp bin get the fraction of the gti duration
                        # that the kp value was within that KP bin

                        # ik selects only those kpwt bins which are non-zero
                        iks = np.where(sawt > 0.0)[0]
                        betabselsa_arr = []
                        bexpotot_arr = []
                        for inds in iks:

                            # FIRST GET THE SOURCE OBSERVATION EXPOSURE in this kp, cor, SA bin

                            # bin ranges for this bin
                            salo = sabin[inds]
                            sahi = sabin[inds + 1]

                            # since the mkf data are sampled once per second, the total time is just the
                            # length of the isel array * 1 second.  Then we have to multiply that time by the
                            # fraction of time the kp value was in that bin, i.e. by kpwt[i]
                            isel = np.where((mkfsasel >= salo) & (mkfsasel <= sahi))[0]

                            # for now let's assume the times are continuous; TODO:check for discontinuities
                            # include Deadtime correction (second term) and fraction of total time in kp,cor bin (third term)
                            sexpo = ((isel * 0.0 + 1.0) * (1.0 - mkfdtsel[isel]) * sawt[inds]).sum()
                            # print(sexpo)
                            if sexpo > gtidur:
                                sexpo = gtidur
                            if verbose:
                                print("    SA: expo = {ex}; dur of GTI = {gtid}".format(ex=sexpo, gtid=gtidur))

                            # SECOND, EXTRACT THE BKG EVENTS and  EXPOSURE in this kp, cor, SA bin

                            # so we've extracted the times when the KPs are within the first bin range;
                            # and the cor's are in their first bin
                            # and are marching through the non-zero sun angle bins
                            # we need to extract BKG events from the bkg events that have KPs and CORs and
                            # SUN_ANGLES in the specified range,
                            # then accumulate the events and the exposure time -
                            # (assume a single "GTI" of duration equal to the exposure time)
                            # betabsel is the background event table in the appropriate KP and COR range
                            betabselsa = betabselc[(betabselc['SUN_ANGLE'] > salo) & (
                                        betabselc['SUN_ANGLE'] <= sahi)]  # number of bkg events in the specified SA bin
                            if verbose:
                                print("    for {0}<Sun Angle<{1}; Number of Events={2}".format(salo, sahi, len(betabselsa)))
                            if len(betabselsa) > 0:
                                betabselsa_arr.append(betabselsa)
                                bhist = np.histogram(betabselc['SUN_ANGLE'], bins=sabin)[0]
                                bexpowt = bhist[inds] / float(
                                    bhist.sum())  # fraction of bkg exposure in this kp, cor, sa bin
                                bdtave = betabselsa[
                                    'MPU_DEADTIME'].mean()  # use the average deadtime for bkg exposure calculation
                                bexpototsa = bexpototc * (1.0 - bdtave) * bexpowt
                                bexpotot_arr.append(bexpototsa)
                            else:
                                print("    No Events Found")
                                bexpototsa = 0.0
                            if verbose:
                                print("    Done SUN Angle index {0}; {1} < SUN_ANGLE < {2}; # Events = {3}".format(inds, salo,
                                                                                                               sahi, len(betabselsa)))
                                print("    Source exposure = {0:.3f}; Background exposure = {1:.3f}".format(sexpo, bexpototsa))
                        btotexpo_kcsa = np.asarray(bexpotot_arr).sum()
                        if verbose:
                            print("    total exposure in this kp, Cor bin for observed SA values = {0}".format(btotexpo_kcsa))
            else:  # gti duration less than minimum
                print("    GTI #{gn} < {gmin} seconds, skipping".format(gn=gtinum, gmin=gtimin))
    # stack all the bkg events in the list of bkg tables
    #print("Length of Events Table = {lt}".format(lt=len(betabselsa_arr)))
    try:
        bkgevt_tot = vstack(betabselsa_arr, join_type='outer')
    except:
        try:
            bkgevt_tot = vstack(betabselsa_arr)
        except:
            print("Could not stack Background Events Array")
            return 0, 0, 0
    bkg_chan = np.arange(1502)
    bkgspec_tot = np.histogram(bkgevt_tot['PI'], bins=bkg_chan)[0]
    # normalize to 50 FPMs
    bkgspec_tot = bkgspec_tot*numfpms/50

    #  Now write an xspec pha file of background
    pha=fits.open(srcpha)
    pha['SPECTRUM'].data['COUNTS'] = np.ceil(bkgspec_tot).astype(int)
    pha['SPECTRUM'].header['EXPOSURE'] = btotexpo_kcsa
    phaout = srcpha.replace('.pha','_bkg.pha')
    if verbose:
        print("\n Writing {po}".format(po=phaout))
    try:
        pha.writeto(phaout, output_verify='fix', checksum=True, overwrite=clobber)
    except Exception as errmsg:
        print("Could not write {out}".format(out=phaout))
        print(errmsg)
    #
    print("Done")
    return bkg_chan[:-1], bkgspec_tot, btotexpo_kcsa


def pfilt_bkgevt(bevtfile, parname='KP', prange=[0,1], clobber=True, verbose=False):
    """
    for an input "enhanced" background event file, create a list of events based on
    a specified range of a specified parameter

    :param bevtfile: name of background event file (generally an UFA file)
    :param parname: name of parameter to filter on
    :param prange: parameter range to filter: accepted events have  prange[0]<= parname < prange[1]
    :param clobber: overwrite if True
    :param verbose: if true write out helpful diagnostics
    :return: astropy table of accepted events
    """
    if not isinstance(bevtfile, (list, tuple, np.ndarray)): # check if bevtfile is a scalar; if so make it a list
        bevtfile=[bevtfile]
    expotot = []
    evtfilt_arr = []
    pmin = prange[0]
    pmax = prange[1]
    for bev in bevtfile:
        if verbose:
            print ("Processing {0}".format(bev))
        evtab = Table.read(bev, hdu='EVENTS', format='fits')
        gti = Table.read(bev, hdu='GTI', format='fits')
        for g in gti:
            expo = 0.0
            # filter by gti first
            gstart = g['START']
            gstop = g['STOP']
            gdur = gstop - gstart
            evtf = evtab[(evtab['TIME'] >= gstart) & (evtab['TIME'] < gstop)]
            totevents = len(evtf)
            if totevents > 0:
                # if gtifilter returned events, then filter by parameter
                evtfilt = evtf[(evtf[parname]>= pmin) & (evtf[parname]<pmax)]
                expo=len(evtfilt)/totevents*gdur
                print('{prange}, {parname} between {pmin:.3f} and {pmax:.3f} Exposure = {expo}'.format(prange=prange, parname=parname, pmin=evtf[parname].min(), pmax=evtf[parname].max(), expo=expo) )
                expotot.append(expo)
                evtfilt_arr.append(evtf)
    evtpfilt = vstack(evtfilt_arr)
    expotot = np.sum(expotot)
    evtpfilt.meta['EXPOSURE'] = expotot
    return evtpfilt

def unit_test(pha='test.pha', obsid=1200040103,root='testdata',
              kpfile = "https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/kp.fits",
              numfpms = 52,
              bevt="https://heasarc.gsfc.nasa.gov/FTP/caldb/data/nicer/xti/pcf/30nov18targskc_enhanced.evt"):
    """
    The unit test will run a sample background calculation.  It assumes the data are located in a subdirectory called "testdata"
    of the current working directory
    :param pha: the name of the pha file
    :param obsid: the test nicer observation id for the dataset
    :param root: the path to the testdata directory
    :param kpfile: The file containing the space-weather KP parameters
    :param bevt: the events file created from the NICER background (blank-sky) observations.
    :return:
    """
    srcpha = os.path.join(root,pha)
    # create mkf3 file from mkf2 file
    mkf2  = os.path.join(root,str(obsid),'auxil','ni1200040103.mkf2')
    mkf3 = mkf2.replace('.mkf2','.mkf3')
    print('Making mkf3 file {0}'.format(mkf3))
    s = add_kp(mkf2, clobber=True, kpfile=kpfile)
    #bkg_chan, bkgspec_tot, bexpotot = mk_bkg_spec_evt(srcpha, gti=gtip, mkffile=mkf3, gtimin=0.1)
    print('Making Background Spectrum')
    bkg_chan, bkgspec_tot, bexpotot = mk_bkg_spec_evt(srcpha, mkf3, bevt=bevt, numfpms = numfpms, gtimin=0, gti=None, clobber=True, verbose=False)
    return bkg_chan, bkgspec_tot, bexpotot

if __name__ == '__main__':
    unit_test()

