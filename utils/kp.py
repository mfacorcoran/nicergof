__author__  = 'Mike Corcoran'
__version__ = '0.3'
__date__    = '2020 Apr 24'



"""
Change History
    2020 Jan 21 Version 0.2 - added update_kpp_def,  mk_kpp_fits
    2020 Apr 24 Version 0.3 - added CVE*, GEOMCOL header keywords

Hey, I have some questions about the Kp file you have in "CALDB," though not actually in CALDB.  I started to Facebook message you about it.

I want to incorporate a Kp value into the prefilter calculator. It would enable alternate COR calculations. I found Nymmik et al 2009, Cosmic Research 47, p191, which has a Kp-modified COR, which maybe useful for our screening.

But... your Kp file that currently "resides" in general CALDB is kind of weird.  It's NICER-specific and not mission generic. (Not indexed but that's not my business).  Also there are alternate sources of Kp, such as Potsdam, which has higher resolution.  I have scripts that produce such a Kp file.

So the upshot is we should agree on a common format so that prefilter can use it.

What I want.
   TIME = MJD (not NICER time)
   KP = KP value (floating point, not integer)
everything else is optional, and not used/needed by prefilter.  A boundary keyword that specifies the origin of the data should be specified.
   CBD10001 = 'ORIGIN("NOAO")'
   CBD10001 = 'ORIGIN("POTSDAM")'

Right now you have a TIME column but it's NICER MET which is kind of mission-specific, and not a appropriate for prefilter, so that's why I ask for alternate time.   It seems we both agree to translate the TIME point to the midpoint between samples, so OK.

Also, you include "-1" values which I think are really "not measured."  Rather than include them I think they should be screened.

Separately I would like to include a separate file for F10.7 for solar activity.  That is also available on-line and needs to be massaged into useful format for prefilter.

Craig

potsdam KP index files can be retrieved from https://www.gfz-potsdam.de/en/kp-index/
i.e. ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/

"""

from ftplib import FTP
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
import os
import subprocess
from subprocess import PIPE
import glob
import numpy as np
from tqdm import tqdm
from nicergof.utils.nicertimeconv import nicertimeconv
from datetime import datetime
import requests

def mk_kp_fits(outfile='kp_noao.fits', header=None,
               DGDdir='/Users/mcorcora/software/github/nicergof/utils/tmp/DGDfiles',
               getdata=True, verbose=False, clobber=True,
               noaoftp="ftp.swpc.noaa.gov", noaoftpdir="/pub/indices/old_indices",
               MJDREFI = 56658, MJDREFF = 0.000777592592592593,LEAPINIT = 2, TIMEZERO = -1):
    """
    mk_kp_fits will create a FITS file containing a history of the interplanetary KP values from NOAO

    NOTE: Uses the NICER ftool nicertimeconv to convert to MET so nicertimeconv needs to be installed

    the kp.fits file is created from the DGD.txt files stored on ftp://ftp.swpc.noaa.gov/pub/indices/old_indices

    The DGD.txt files are of the form:
    #                Middle Latitude        High Latitude            Estimated
    #              - Fredericksburg -     ---- College ----      --- Planetary ---
    #  Date        A     K-indices        A     K-indices        A     K-indices
    2017 01 01     9  3 3 3 2 1 2 1 1    26  2 4 4 6 5 2 2 1    14  4 4 3 2 2 3 2 2

    The columns needed are the Estimated Planetary K indices.
    Each day has 8 of them since they are reported every 3 hours.
    KCG uses the following 8 time assignments:
    T01:30:00, T04:30:00,T07:30:00,T10:30:00,T13:30:00,T16:30:00,T19:30:00,T22:30:00

    :parameter outfile: name of output kp.fits file (with directory path if needed)
    :parameter header: template FITS header for kp table
    :parameter DGDdir: directory which holds the DGD files from noao
    :parameter getdata: if True, will retrieve data
    :parameter verbose: if True, will print diagnostic messages
    :parameter clobber: if True, will overwrite existing outfile
    :parameter noaoftp: NOAO FTP site
    :parameter noaoftpdir: NOAO directory in FTP site where the DGD files are located
    :return: returns a status (= 0 if no errors; otherwise a diagnostic message)

    Changes:
    MFC 20200115 ADDED CBD10001 = 'NOAO' (request from C. Markwardt) & TNULL2 = -1 to  denote null values;
            also added TIMEPIXR keyword and CVSD, CVST values

    """
    # assume everything will work
    status = 0
    skipit = False

    # if getdata is True then download the data from the NOAO ftp site to the DGDdir
    if getdata:
        if verbose:
            print('Getting KP data from NOAA')
        try:
            get_kp_data(DGDdir=DGDdir, noaoftp=noaoftp, noaoftpdir=noaoftpdir, verbose=verbose)
        except Exception as errmsg:
            status = "Error getting data from ftp://{0}{1} ({2})".format(noaoftp, noaoftpdir, errmsg)
            print(status)
            return status
    # read the data from each DGD file
    ## setup some lists and variable definitions

    MJD = []
    KP = []
    kptime = ["T01:30:00",
              "T04:30:00",
              "T07:30:00",
              "T10:30:00",
              "T13:30:00",
              "T16:30:00",
              "T19:30:00",
              "T22:30:00"]
    counter = range(0, 17, 2)

    ## get list of dgd files
    dgd = glob.glob(os.path.join(DGDdir, '*DGD*.txt'))
    if len(dgd) == 0:
        status = "No DGD files in {0}".format(DGDdir)
        skipit = True
    # just get data since year 2000, i.e. files that start with 2*_DGD.txt
    dgd = [x for x in dgd if '2' in os.path.split(x)[-1]]
    if not skipit:
        for fname in dgd:
            #fname = os.path.join(DGDdir, d)
            print("File = {0}".format(fname))
            with open(fname, 'r') as f:
                for line in f:
                    if '20' in line[0:2]:
                        date = []
                        kpval = []
                        # parse line
                        kpl = line[63:]
                        for i in counter[:-1]:
                            kpv = kpl[i:i + 2].strip()
                            try:
                                kpval.append(int(kpv))
                            except Exception as e:
                                print(fname, e)
                        newdate = line[0:10].replace(' ', '-')
                        for kt in kptime:
                            t = "{date}{kt}".format(date=newdate, kt=kt)
                            date.append(Time(t, format='fits').mjd)
                        MJD.extend(date)
                        KP.extend(kpval)
        # sort by MJD
        iso = np.argsort(MJD)
        MJD = np.asarray(MJD)[iso]
        KP = np.asarray(KP)[iso]
        UTC = [Time(x,format='mjd').isot for x in MJD]
        print("starting met calculation")
        met = nicertimeconv(MJD,informat='mjd', outformat='met')
        print("... finished MET calculation")
        c1 = fits.Column(name='TIME', array=MJD, format='1D')
        c2 = fits.Column(name='KP', array=KP, format='J')
        c3 = fits.Column(name='MET', array=met, format='1D')
        c4 = fits.Column(name='UTC', array=UTC, format='23A')
        # update header with CALDB keywords
        if header == None:
            header = fits.Header()
        header['EXTNAME'] = 'NOAO_KP'
        header['ORIGIN'] = "NOAO"
        header['HDUCLAS'] = 'OGIP'
        header['HDUCLAS1'] = 'TEMPORALDATA'
        header['CREATOR'] = 'kp.py'
        tnow = Time.now().fits
        header['DATE'] = tnow
        header['MJDREFI'] = MJDREFI
        header.comments['MJDREFI'] = '2014-01-01T00:00:00'
        header['MJDREFF'] = MJDREFF
        header['LEAPINIT'] = LEAPINIT
        header.comments['LEAPINIT'] = 'Leap seconds since MJDREFI'
        header['TIMEZERO'] = TIMEZERO
        header['TNULL2'] = -1
        header.comments['TNULL2'] = 'KP not measured'
        header['TIMEPIXR'] = 0.5
        header.comments['TIMEPIXR'] = 'Times refer to the center of the time bin'
        header['GEOMCOL'] = 'KP'
        header.comments['GEOMCOL'] = 'Name of primary geomagnetic column'
        cvsd = Time(MJD.min(), format='mjd').iso[0:10]
        cvst = Time(MJD.min(), format='mjd').iso[11:].split('.')[0]
        cved = Time(MJD.max(), format='mjd').iso[0:10]
        cvet = Time(MJD.max(), format='mjd').iso[11:].split('.')[0]
        calkeywords = {'TELESCOP': 'GEN', 'INSTRUME': 'INS', 'DETNAM': '', 'FILTER': '',
                   'CCLS0001': 'PCF', 'CCNM0001': 'KP_VALUES', 'CDTP0001': 'DATA',
                   'CVSD0001': '{0}'.format(cvsd), 'CVST0001': '{0}'.format(cvst), 'CBD10001': 'SOURCE("NOAO")',
                   'CBD20001': '', 'CBD30001': '', 'CBD40001': '', 'CBD50001': '',
                   'CBD60001': '', 'CBD70001': '', 'CBD80001': '', 'CBD90001': '',
                   'CDES0001': 'Interplanetary KP values from NOAO'}
        for k in calkeywords.keys():
            header[k] = calkeywords[k]
        header['CVED0001'] = cved
        header.comments['CVED0001'] = 'Validity end date (non-standard CALDB keyword)'
        header['CVET0001'] = cvet
        header.comments['CVET0001'] = 'Validity end time (non-standard CALDB keyword)'
        header.add_comment('KP data from ftp://ftp.swpc.noaa.gov/pub/indices/old_indices')
        header.add_comment('FITS file created by mk_kp_fits()')
        header.add_comment('On {0} UT'.format(tnow))
        t = fits.BinTableHDU.from_columns([c1, c2, c3, c4], header=header)
        t.header.comments['TTYPE1'] = 'Time in MJD'
        t.header.comments['TTYPE2'] = 'KP value'
        t.header.comments['TTYPE3'] = 'Time of KP measurement in NICER MET'
        t.header.comments['TTYPE4'] = 'Time in ISOT format'
        try:
            t.writeto(outfile, output_verify='exception', overwrite=clobber, checksum=True)
        except Exception as e:
            status = "Error writing {0} ({1})".format(outfile, e)
            print(e)
    if status != 0:
        print(status)
    else:
        print('Finished writing {0}'.format(outfile))
    return status


def get_kp_data(DGDdir='/software/github/nicergof/kcg/Notebooks/resources/DGDfiles',
                noaoftp="ftp.swpc.noaa.gov", noaoftpdir="/pub/indices/old_indices",
                verbose=False):
    """
    Gets the DGD files from the NOAO ftp site and stores them locally
    :parameter DGDdir: directory containing the downloaded DGD files
    :parameter noaoftp: the noao ftp site
    :parameter noaoftpdir: the directory on the noao site which holds the DGD files
    :parameter verbose: if True prints diagnostic messages
    :return: no return value

    updates:
      mfc 20190214: fixed bug if DGD files not previously downloaded (mike l report)
    """
    # get list of files and download them using wget
    dgddir = DGDdir

    # Get file list
    ftp = FTP(noaoftp)
    ftp.login()
    ftp.cwd(noaoftpdir)

    data = []
    ftp.dir(data.append)
    ftp.close()

    dgd = os.listdir(DGDdir)
    if len(dgd) == 0:
        dgd=[]

    for d in data:
        #if ('DGD' in d) and ('Q' in d):
        if ('DGD' in d):
            test = d.split()[-1]
            #if int(test[0:4]) >= 2017:
            #    dgd.append(test)
            dgd.append(test)

    for d in dgd:
        fname = os.path.join(dgddir, d)
        # if file does not exist, get it
        if not os.path.isfile(fname):
            cmd = 'wget --output-document={0} ftp://{1}/{2}/{3}'.format(fname, noaoftp, noaoftpdir, d)
            if verbose:
                print(cmd)
            os.system(cmd)
        else:
            if verbose:
                print("{0} exists".format(fname))

    # update the most current file
    curfile = dgd[-1]
    print("Updating {0}".format(curfile))
    fname = os.path.join(dgddir, curfile)
    cmd = 'wget --output-document={0} ftp://{1}/{2}/{3}'.format(fname, noaoftp, noaoftpdir, d)
    if verbose:
        print(cmd)
    os.system(cmd)
    return

def get_kpp_tab(kppftp = "ftp.gfz-potsdam.de", kppdir = "pub/home/obs/kp-ap/tab",
                kptabdir = '/Users/mcorcora/software/github/nicergof/utils/refdata/potsdam',
                verbose=True, get_all=False):
    """
    Downloads the Potsdam KP files from the kp*.tab files on the Potsdam KP FTP site
    ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/
    (also see https://www.gfz-potsdam.de/en/kp-index/ for more information)

    :parameter kppdef: output name of the "definitive" kp file
    :parameter kppftp: ftp site
    :parameter kppdir: directory where the kp*tab files are located
    :return:
    """
    status = 0

    # Get file list
    ftp = FTP(kppftp)
    ftp.login()
    ftp.cwd(kppdir)

    data = []
    ftp.dir(data.append)
    ftp.close()

    kpfile = [x.split()[-1] for x in data if '.tab' in x]
    for k in kpfile:
        fname = os.path.join(kptabdir,k)
        if get_all:
            cmd = 'wget --output-document={0} ftp://{1}/{2}/{3}'.format(fname, kppftp, kppdir, k)
            if verbose:
                print(cmd)
            os.system(cmd)
        else:
            if not os.path.isfile(fname):
                cmd = 'wget --output-document={0} ftp://{1}/{2}/{3}'.format(fname, kppftp, kppdir, k)
                if verbose:
                    print(cmd)
                os.system(cmd)
            else:
                if verbose:
                    print("{0} already exists".format(fname))
    return status

def mk_kpp_deffile(kppdeffile='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_def_2019.txt',
                   kpptabdir = '/software/github/nicergof/utils/refdata/potsdam',
                   dformat='%y%m%d'):
    """
    From kp.tab files in kptabdir, create the "definitive" kp ascii file
    :parameter kppdeffile: the "definitive KP file" to create (default /software/github/nicergof/utils/refdata/potsdam/mfc/kpp_definitive.txt)
    :parameter kptabdir: directory where the downloaded Potsdam KP files are stored
    :return:
    """
    status = 0
    print('Creating {0}'.format(kppdeffile))
    kppfiles = os.listdir(kpptabdir)
    kppfiles = [k for k in kppfiles if '.tab' in k]
    kppdat = []
    for k in kppfiles:
        with open(os.path.join(kpptabdir, k), 'r') as f:
            ll = f.readlines()
            dat = [x[:32] for x in ll if len(x[4:5].strip()) != 0]
        kppdat.extend(dat)
    tkppdat = [Time(datetime.strptime(x[:6], dformat)).mjd for x in kppdat]
    isort = np.argsort(tkppdat)
    # sort the kp data array by time (mjd)
    defkpp = np.asarray(kppdat)[isort]
    with open( kppdeffile, 'w') as f:
        for l in defkpp:
            f.write(l+'\n')
    return status

def kpp_read_def(kppdef = '/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_definitive.txt',
                 dformat='%y%m%d'):
    print('Reading Definitive KPP file {0}'.format(kppdef))
    kptbinmid = np.arange(0, 24, 3) + 1.5
    kpdict = {'o': 0.0, '+': 1. / 3., '-': -1. / 3.}
    kp=[]
    kpmjd=[]
    kp_orig=[]
    with open(kppdef, 'r') as f:
        ll = f.readlines()
        for l in ll:
            kparr = l[:32].split()
            kpval = kparr[1:]
            kpdt = datetime.strptime(kparr[0], dformat)
            for t, k in zip(kptbinmid, kpval):
                kp.append(int(k[0]) + kpdict[k[1]])
                kpmjd.append(Time(kpdt).mjd + t / 24.)
                kp_orig.append(k)
        kpdef = {'KP':kp, 'KP_orig':kp_orig, 'MJD':kpmjd}
        return kpdef


def update_kpp_def(kppout='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_current.txt',
        kppdef = '/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_def_2019.txt',
        pkpql = "http://www-app3.gfz-potsdam.de/kp_index/pqlyymm.tab",
        dformat = '%y%m%d', header=None):
    """
    This routine reads the definitive KP file (as created by mk_kpp_deffile), reads the latest KPP data file, and
    appends the two, and outputs to the file given by the kppout parameter.

    The Potsdam files used to create the "definitive file" are given in the kp*.tab files in
    ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/

    The columns in the file are
    MJD, ISO time, KP values from the Potsdam data tables, KP value converted to floats

    :parameter kppout: name of the output file which contains the combined "definitive" data with the current month's data appended.
    :parameter kpdef: name of the ascii file containing the definitive list
    :parameter pkpql: url giving the quick-look data for the current month
    :return: status = 0 if completion without error
    """
    status=0
    kptbinmid = np.arange(0, 24, 3) + 1.5
    kpdict = {'o': 0.0, '+': 1. / 3., '-': -1. / 3.}
    print('Reading Definitive KPP file {0}'.format(kppdef))
    with open(kppdef, 'r') as f:
        ll = f.readlines()
    with open(kppout, 'w') as f:
        print('Updating {0}'.format(kppout))
        for l in ll:
            kparr = l[:32].split()
            kpval = kparr[1:]
            kpdt = datetime.strptime(kparr[0], dformat)
            for t, k in zip(kptbinmid, kpval):
                kp= (int(k[0]) + kpdict[k[1]])
                kpmjd = (Time(kpdt).mjd + t / 24.)
                kp_orig= k
                # kp, kpval, dt, Time(dt,format='mjd').iso
                outstring = "{3:.4f} {0} {1} {2:.4f}".format(Time(kpmjd, format='mjd').iso,k,kp,kpmjd)
                #print(outstring)
                f.write(outstring + "\n")
        # add QL data
        print('Adding data for current month to {0}'.format(kppout))
        req = requests.get(pkpql)
        pkpqldata = req.text.split('\n')
        pkpqldata = [x for x in pkpqldata if len(x)>0] # get rid of blank lines
        for l in pkpqldata:
            kparr = l[:32].split()
            kpval = kparr[1:]
            try:
                kpdt = datetime.strptime(kparr[0], dformat)
            except Exception as errmsg:
                print(errmsg)
                status = -1
                return status
            for t, k in zip(kptbinmid, kpval):
                skipit=False
                try:
                    kp = int(k[0]) + kpdict[k[1]]
                except Exception as errmsg:
                    print('Error writing QL data: ({0})', format(errmsg))
                    skipit=True
                if not skipit:
                    kpmjd = (Time(kpdt).mjd + t / 24.)
                    outstring = "{0:.4f} {1} {2} {3:.4f}".format(kpmjd,Time(kpmjd, format='mjd').iso,k,kp)
                    f.write(outstring + "\n")
    return status

def mk_kpp_fits(kppcurfile='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_current.txt',
                kppfits = '/software/github/nicergof/utils/refdata/potsdam/mfc/kpp.fits',
                overwrite=True):
    """
    This function writes the potsdam kp data in fits format

    :param kppcurfile: name of the ascii file containing the historical and latest potsdam kp data
    :param kppfits: name of the fits output file containing the latest potsdam kp data
    :return: status (0 if no error)
    """
    status = 0
    print('Reading {0}'.format(kppcurfile))
    kptab = Table.read(kppcurfile,
                     names=['MJD', 'DATE', 'TIME', 'KP_ORIG', 'KP'], format='ascii')
    kptab.meta['EXTNAME'] = 'KP_POTSDAM'
    kptab.meta['HDUCLAS'] = 'OGIP'
    kptab.meta['HDUCLAS1'] = 'TEMPORALDATA'
    kptab.meta['CREATOR'] = 'kpp from kp.py'
    tnow = Time.now().fits
    kptab.meta['DATE'] = tnow
    cvsd = kptab['DATE'][0]
    cvst = kptab['TIME'][0]
    calkeys = {'TELESCOP': 'GEN', 'INSTRUME': 'INS', 'DETNAM': '', 'FILTER': '',
               'CCLS0001': 'PCF', 'CCNM0001': 'KP_VALUES', 'CDTP0001': 'DATA',
               'CVSD0001': '{0}'.format(cvsd), 'CVST0001': '{0}'.format(cvst), 'CBD10001': 'SOURCE("Potsdam")',
               'CBD20001': '', 'CBD30001': '', 'CBD40001': '', 'CBD50001': '',
               'CBD60001': '', 'CBD70001': '', 'CBD80001': '', 'CBD90001': '',
               'CDES0001': 'Interplanetary KP values from Potsdam'}
    for k in calkeys:
        kptab.meta[k] = calkeys[k]
    kptab.meta['TIMEPIXR'] = 0.5
    print('Writing {0}'.format(kppfits))
    try:
        kptab.write(kppfits,format='fits',overwrite=overwrite)
    except:
        print("Error writing output file {0}".format(kppfits))
        status = -1
    return status


if __name__ == '__main__':
    outfile = '/Users/mcorcora/software/github/nicergof/utils/tmp/kp_noao.fits'
    caldir = "/Users/mcorcora/software/github/nicergof/utils/tmp"
    dgddir = "/Users/mcorcora/software/github/nicergof/utils/tmp/DGDfiles"
    ftpsite = "ftp.swpc.noaa.gov"
    rdir = "/pub/indices/old_indices"
    # mk_kp_fits(outfile=outfile, header=None, DGDdir=dgddir,
    #          getdata=False, verbose=True, clobber=True,
    #           noaoftp=ftpsite, noaoftpdir=rdir)
    # print("starting")
    mk_kp_fits(outfile=outfile, DGDdir=dgddir)
    # update_kpp_def(kppout='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_current.txt',
    #     kppdef = '/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_def_2019.txt',
    #     pkpql = "http://www-app3.gfz-potsdam.de/kp_index/qlyymm.tab",
    #     dformat = '%y%m%d', header=None)

