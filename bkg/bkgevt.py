import numpy as np
import pandas as pd
from subprocess import Popen, PIPE
from bs4 import BeautifulSoup
import os
from astropy.table import Table, vstack, Column
from astropy.io import fits
from astropy.time import Time
from nicergof.bkg.bkg_estimator import add_kp
import glob

"""
    METHOD:
    a) get list of event files from bkg observations 
    b) get the mkf files for each event file
    c) read the events file as an astropy table
    d) get the KP, COR_SAX, SUN_ANGLE and MPU_DEADTIME columns from the MKF file
    e) add the KP, COR_SAX, SUN_ANGLE and MPU_DEADTIME columns (interpolated to the event times) to the events file table
    f) merge the event file tables using vstack to join the tables by rows
    g) sort the merged table by time

"""

def find_bkg_data(obsfile = 'BKGD_RXTE_obsids.txt', nicer_archive='/FTP/nicer/data/obs', datatype='mpu7_ufa'):
    """
    scans the NICER archive and generates a list of background event files
    for the BKGD_RXTEn fields
    
    :parameter obsfile: ascii file that contains the list of unique obsid identifiers for the BKGD_RXTEn fields;
    should be pipe-delimited and have columns which include the Field name and OBSID, with no header line
    :parameter nicer_archive: location of the NICER data archive
    """
    # get list of arhive directories
    #dirs = Popen('ls -d {0}/20*'.format(nicer_archive), stdout=PIPE, shell=True, universal_newlines=True).communicate()[0].split()
    # # get bkg root obsids
    # with open(obsfile,'r') as f:
    #     obsroots = f.readlines()
    # obsroots = [x.strip().split()[1] for x in obsroots]
    bkgobs = []
    df = pd.read_csv(obsfile, sep='|', usecols=[1, 2], names=['Field', 'OBSID'])
    obsroots =["{:010d}".format(x) for x in df.OBSID]
    for o in obsroots:
        cmd = 'ls -d  /FTP/nicer/data/obs/20*/{0}*'.format(o)
        b = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True).communicate()[0].split()
        bkgobs.extend(b)
    evtlist=[]
    for b in bkgobs:
        ef = os.path.join(b,'xti/event_cl','*{0}*'.format(datatype))
        cmd = 'ls -d {0}'.format(ef)
        ev = Popen(cmd,shell=True, stdout=PIPE, universal_newlines=True).communicate()[0].split()
        evtlist.extend(ev)
    return evtlist
    
    
def mk_evtbk_table(bevtfile,outfile='compose', indir=None,scratchdir='/tmp',
                   kp='/FTP/caldb/data/gen/pcf/kp.fits',
                   mkcols=['KP','SUN_ANGLE','COR_SAX','MPU_DEADTIME'], clobber=True,
                   memmap=True):
    """
    For a given bkg event file, find the mkf2 file (or create it if doesn't exist) for the event file, 
    create a table with additional information from the mkf file (including KP) 
    then stack the event tables and sort them by time
    
    :parameter bevtfile: FITS background event file
    :parameter outfile: name of output "enhanced" background FITS file with additional info from mkf3 file included; if None, don't write output file; if "compose" create background file with standard name in scratchdir
    :parameter indir: location of background event files; if None, use the NICER archive
    :parameter scratchdir: location of temporary and output file
    :parameter kp: name of the FITS file containing the KP solar activity data
    :return: combined events table
    
    """
    mkf =  bevtfile.replace('xti/event_cl','auxil').replace('_0mpu7_ufa.evt.gz','.mkf.gz')
    mkffile = os.path.split(mkf)[-1]
    if not indir:
        indir = os.path.split(os.path.split(mkf)[0])[0]
        print('INDIR = {0}'.format(indir))
    # make mkf2 file in scratch area then add kp
    mkf2 = os.path.join(scratchdir, mkffile.replace('.mkf','.mkf2').replace('.gz',''))
    mkf3 = mkf2.replace('mkf2','mkf3')
    if not os.path.isfile(mkf3):
        print(mkf2)
        cmd = 'niprefilter2 indir = {0} infile={1} outfile = {2} clobber=yes'.format(indir, mkf, mkf2)
        print(cmd)
        os.system(cmd)
        # add KP to mkf2 file
        stat=add_kp(mkf2,kpfile=kp,clobber=True)
    # read mkf3 file as Table
    print('Reading {0}'.format(mkf3))
    mkf3dat = fits.open(mkf3)['PREFILTER'].data
    # read bkg evt file into an astropy Table & sort by time
    print('Reading {0}'.format(bevtfile))
    evtab = Table.read(bevtfile,format='fits', hdu='EVENTS', memmap=memmap)
    evtab.sort(['TIME'])
    # add OBSID as a column in the evtab
    obsid = os.path.split(mkf3)[-1][2:12]
    colobs = Column(name='OBS_ID',data=[obsid]*len(evtab))
    evtab.add_column(colobs)
    # get times, Sun angle, kp, cor_sax and deadtime from mkf3 file
    skip=False
    for m in mkcols:
        if len(mkf3dat[m].shape) == 2:
            mkvals = mkf3dat[m].mean(axis=1)
        else:
            mkvals = mkf3dat[m]
        try:
            coli = np.interp(evtab['TIME'],mkf3dat['TIME'],mkvals)
        except Exception as errmsg:
            print("Error ({0}) getting column {1} from {0}".format(errmsg, m, mkf3))
            skip=True
        if not skip:
            newcol = Column(name=m,data=coli)
            evtab.add_column(newcol)
    # convert COR_SAX, KP, SUN_ANGLE, MPU_DEADTIME from 64bit to 32 bits
    f32cols = ['COR_SAX', 'KP', 'SUN_ANGLE','MPU_DEADTIME']
    for f in f32cols:
        evtab[f]=np.float32(evtab[f])
    if outfile:
        if outfile == 'compose':
            outfile = os.path.split(mkf3)[-1].replace('.mkf3','_enhanced.evt')
            outfile = os.path.join(scratchdir, outfile)
        gtitab = Table.read(bevtfile, hdu='GTI', format='fits')
        tablist = [evtab, gtitab]
        print('Writing {0}'.format(outfile))
        #evtab.write(outfile, format='fits', overwrite=clobber)
        mk_fits_table(tablist,outfile, clobber=clobber)
    return evtab, gtitab
    
    
def combine_bevts(bevtlist, outfile='compose', scratchdir='/tmp'):
    """
    For a list of "enhanced" background event files combine events files

    Issues:
        - need to have the merged event file sorted in time
        - need to update exposure time in the combined table
        - need to create a GTI extension which holds the start/stop time for each BKG field obsid (columns
        OBSID, START, STOP, DURATION)
    Method
        - use vstack to merge the enhanced bkg event lists for each obsid 
        - update the merged table .meta['EXPOSURE'] for the total exposure time
        - update the GTI extension for the start/stop times of each obsid
        - create fits.Column objects from the columns in the merged table
        - create an hdu from the Column objects using fits.BinTableHDU.from_columns()
        - create an HDULIST with the merged event tables and the good time intervals, then write out
        
    :parameter evtlist: a list of the names of the enhanced bkg evt files
    :parameter outfile: name of the output file (with path)
    
    :return: a 0 status if all went ok, non-zero otherwise
    
    """
    stat=0
    print('Reading {0}'.format(bevtlist[0]))
    evtot = Table.read(bevtlist[0],format='fits', hdu='EVENTS',  memmap=True)
    print('Number of events ={0}'.format(len(evtot)))
    gtitot = Table.read(bevtlist[0],format='fits',hdu='GTI')
    gtitot['Duration']=gtitot['STOP']-gtitot['START']
    mjdref = gtitot.meta['MJDREFF'] + gtitot.meta['MJDREFI']
    gtitot['STARTISO'] = Time(gtitot['START'].data / 86400 +mjdref, format='mjd').iso
    gtitot['STOPISO'] = Time(gtitot['STOP'].data / 86400 + mjdref, format='mjd').iso
    gtitot['EVT_FILE'] = os.path.split(bevtlist[0])[-1]
    expotot=np.sum(gtitot['Duration'])
    # tstart = Time(evtot.meta['DATE-OBS']).mjd
    # tstop  = Time(evtot.meta['DATE-END']).mjd
    # expotot = evtot.meta['EXPOSURE']
    # print('Expo = {0}'.format(expotot))
    # gti = {
    #     'START': [tstart],
    #     'STOP': [tstop],
    #     'Duration':[(tstop - tstart)*86400.0],
    #     'OBSID':[evtot['OBS_ID'][0]],
    #     'EXPOSURE':[expotot]}
    evtot = rem_bevt_keys(evtot)
    evlist=[evtot]
    gtilist=[gtitot]
    for bev in bevtlist[1:]:
        print('Adding {0}'.format(bev))
        evt = Table.read(bev,format='fits', hdu='EVENTS', memmap=True)
        print('Number of events ={0}'.format(len(evt)))
        tstart = Time(evt.meta['DATE-OBS']).mjd
        tstop  = Time(evt.meta['DATE-END']).mjd
        gti = Table.read(bev,format='fits',hdu='GTI')
        gti['Duration'] = gti['STOP'] - gti['START']
        gti['STARTISO'] = Time(gti['START'].data / 86400 + mjdref, format='mjd').iso
        gti['STOPISO'] = Time(gti['STOP'].data / 86400 + mjdref, format='mjd').iso
        gti['EVT_FILE'] = os.path.split(bev)[-1]
        # gti['START'].append(tstart)
        # gti['STOP'].append(tstop)
        # gti['Duration'].append((tstop - tstart)*86400.0)
        # gti['OBSID'].append(evt['OBS_ID'][0])
        # gti['EXPOSURE'].append(evt.meta['EXPOSURE'])
        # expo = evt.meta['EXPOSURE']
        expo = np.sum(gti['Duration'])
        print('Expo = {0}'.format(expo))
        expotot = expotot + expo
        print('Total exposure= {0}'.format(expotot))
        evt=rem_bevt_keys(evt)
        #evlist.append(evt)
        evtot=vstack([evtot, evt],metadata_conflicts='silent')
        del evt # delete evt table to save memory
        gti=rem_bevt_keys(gti)
        #gtilist.append(gti)
        gtitot = vstack([gtitot,gti],metadata_conflicts='silent')
    #evtot=vstack(evlist,metadata_conflicts='silent')
    #gtitot=vstack(gtilist,metadata_conflicts='silent')
    evtot.sort(['TIME']) 
    evtot.meta['EXPOSURE']=expotot
    #gtitab = Table(gti)
    #gtitab.meta['EXTNAME']='GTI'
    #gtitab.sort(['START'])
    gtitot.meta['EXTNAME'] = 'GTI'
    gtitot.sort(['START'])
    if outfile:
        if outfile == 'compose':
            outfile = os.path.join(scratchdir, 'nicer_combined_bkg_enhanced.evt')
        tablist=[evtot, gtitot]
        try:
            mk_fits_table(tablist,outfile, clobber=True)
        except Exception as errmsg:
            print('Problem writing {0} ({1}); returning'.format(outfile, errmsg))
            stat=1
    del evtot # delete to save memory
    return stat
    
def rem_bevt_keys(evtab):
    """
    convenience function to remove useless keys from the bkg event table before combining to avoid annoying warnings
    """
    try:
        del evtab.meta['DATE']
    except:
        pass
    try:
        del evtab.meta['ONTIME']
    except:
        pass
    try:
        del evtab.meta['OBS_ID']
    except:
        pass
    try:
        del evtab.meta['EXPOSURE']
    except:
        pass
    try:
        del evtab.meta['CHECKSUM']
    except:
        pass
    try:
        del evtab.meta['DATASUM']
    except:
        pass
    try:
        del evtab.meta['RA_NOM']
    except:
        pass
    try:
        del evtab.meta['DEC_NOM']
    except:
        pass
    try:
        del evtab.meta['TSTART']
    except:
        pass
    try:
        del evtab.meta['TSTOP']
    except:
        pass
    try:
        del evtab.meta['DATE-OBS']
    except:
        pass
    try:
        del evtab.meta['DATE-END']
    except:
        pass
    try:
        del evtab.meta['TELAPSE']
    except:
        pass
    return evtab
                          
    
def mk_fits_table(tablist,outfile, clobber=True):
    """
    From a list of astropy tables, create a FITS file with the tables as bintable extensions
    :parameter tablist: list of astropy tables
    :parameter outfile: name of output fits file
    """
    stat=0
    hdulist=[fits.PrimaryHDU()]
    for t in tablist:
         hdulist.append(fits.BinTableHDU(t))
    hdu = fits.HDUList(hdulist)
    print('Writing {0}'.format(outfile))
    try:
        hdu.writeto(outfile, overwrite=clobber, output_verify='fix')
    except Exception as errmsg:
        print('Problem writing {0} {1}; returning'.format(outfile, errmsg))
        stat = 1
    return stat

if __name__ == '__main__':
    bevtfile = glob.glob("/Users/mcorcora/SXDC/Data/NICER/kp_model/20181130/work/ni*")
    print(bevtfile)
    for bf in bevtfile:
        evtPF = pfilt_bkgevt(bf, parname='KP', prange=[0,2], clobber=True)
        print(evtPF.meta['EXPOSURE'])