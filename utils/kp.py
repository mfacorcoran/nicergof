

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

def mk_kp_fits(outfile='kp_test.fits', header=None,
               DGDdir='/software/github/nicergof/kcg/Notebooks/resources/DGDfiles',
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

    :param outfile: name of output kp.fits file (with directory path if needed)
    :param header: template FITS header for kp table
    :param DGDdir: directory which holds the DGD files from noao
    :param getdata: if True, will retrieve data
    :param verbose: if True, will print diagnostic messages
    :param clobber: if True, will overwrite existing outfile
    :param noaoftp: NOAO FTP site
    :param noaoftpdir: NOAO directory in FTP site where the DGD files are located
    :return: returns a status (= 0 if no errors; otherwise a diagnostic message)

    """
    # assume everything will work
    status = 0
    skipit = False

    # if getdata is True then download the data from the NOAO ftp site to the DGDdir
    if getdata:
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
                            kpval.append(int(kpv))
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
        #
        # Create kp.fits file - TODO: calculate MET time from MJD
        #
        # met = []
        # ind = arange(0,5).astype(int)
        print("starting met calculation")
        #for i, m in tqdm(enumerate(MJD)):
        # for m in tqdm(MJD):
        #     #t = subprocess.run(['nicertimeconv', '-met', Time(m, format='mjd').isot], universal_newlines=True, stdout=PIPE).stdout
        #     t = subprocess.Popen(['nicertimeconv', '-met', Time(m, format='mjd').isot], universal_newlines=True, stdout=PIPE).communicate()[0]
        #     met.append(float(t))
        #     # print(float(t), m, Time(m,format='mjd').isot)
        # met = np.asarray(met)
        met = nicertimeconv(MJD,informat='mjd', outformat='met')
        print("... finished MET calculation")
        # c1 = fits.Column(name='TIME', array=np.zeros(len(KP)), format='1D')
        c1 = fits.Column(name='TIME', array=met, format='1D')
        c2 = fits.Column(name='KP', array=KP, format='J')
        c3 = fits.Column(name='MJD', array=MJD, format='1D')
        c4 = fits.Column(name='UTC', array=UTC, format='23A')
        # update header with CALDB keywords
        if header == None:
            header = fits.Header()
        header['EXTNAME'] = 'KP'
        header['HDUCLAS'] = 'OGIP'
        header['HDUCLAS1'] = 'TEMPORALDATA'
        header['CREATOR'] = 'kp.py'
        header['MJDREFI'] = MJDREFI
        header.comments['MJDREFI'] = '2014-01-01T00:00:00'
        header['MJDREFF'] = MJDREFF
        header['LEAPINIT'] = LEAPINIT
        header.comments['LEAPINIT'] = 'Leap seconds since MJDREFI'
        header['TIMEZERO'] = TIMEZERO
        calkeys = {'TELESCOP': 'NICER', 'INSTRUME': 'XTI', 'DETNAM': '', 'FILTER': '',
                   'CCLS0001': 'CPF', 'CCNM0001': 'KP_VALUES', 'CDTP0001': 'DATA',
                   'CVSD0001': '2017-01-01', 'CVST0001': '00:00:00', 'CBD10001': '',
                   'CBD20001': '', 'CBD30001': '', 'CBD40001': '', 'CBD50001': '',
                   'CBD60001': '', 'CBD70001': '', 'CBD80001': '', 'CBD90001': '',
                   'CDES0001': 'Interplanetary KP values from NOAO'}
        for k in calkeys:
            header[k] = calkeys[k]
        header.add_comment('KP data from ftp://ftp.swpc.noaa.gov/pub/indices/old_indices')
        header.add_comment('FITS file created by mk_kp_fits()')
        header.add_comment('On {0} UT'.format(Time.now().isot))
        t = fits.BinTableHDU.from_columns([c1, c2, c3, c4], header=header)
        t.header.comments['TTYPE1'] = 'Time of KP measurement in NICER MET'
        t.header.comments['TTYPE2'] = 'KP value'
        t.header.comments['TTYPE3'] = 'Time in MJD'
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
    :param DGDdir: directory containing the downloaded DGD files
    :param noaoftp: the noao ftp site
    :param noaoftpdir: the directory on the noao site which holds the DGD files
    :param verbose: if True prints diagnostic messages
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
        if ('DGD' in d) and ('Q' in d):
            test = d.split()[-1]
            if int(test[0:4]) >= 2017:
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
                kptabdir = '/software/github/nicergof/utils/refdata/potsdam',
                verbose=True, get_all=False):
    """
    Creates the kppdef ascii file from the kp*.tab files on the Potsdam KP FTP site
    ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/

    :param kppdef: output name of the "definitive" kp file
    :param kppftp: ftp site
    :param kppdir: directory where the kp*tab files are located
    :return:
    """

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
    return

def mk_kpp_deffile(kppdeffile='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_definitive.txt',
                   kpptabdir = '/software/github/nicergof/utils/refdata/potsdam',
                   dformat='%y%m%d'):
    """
    From kp.tab files in kptabdir, create "definitive" kp ascii file
    :param kppdeffile:
    :param kptabdir:
    :return:
    """
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
    with open(os.path.join(kpptabdir, kppdeffile),'w') as f:
        for l in defkpp:
            f.write(l+'\n')
    return



def kpp(kppout='/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_current.txt',
        kppdef = '/software/github/nicergof/utils/refdata/potsdam/mfc/kpp_definitive.txt',
        pkpql = "http://www-app3.gfz-potsdam.de/kp_index/qlyymm.tab",
        dformat = '%y%m%d'):
    """
    This routine creates a fits version of the Potsdam KP data.  It includes the "definitive" data given in the kpkcg file
    along with the quicklook data for the current month.

    The definitive files are given in the kp*.tab files in
    ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/

    The columns in the file are
    MJD, ISO time, KP values from the Potsdam data tables, KP value converted to floats

    :param kppout: name of the output file
    :param kpdef: name of the ascii file containing the definitive list
    :param pkpql: url giving the quick-look data for the current month
    :return: status = 0 if completion without error
    """
    print('Writing {0}'.format(kppout))
    kptbinmid = np.arange(0, 24, 3) + 1.5
    kpdict = {'o': 0.0, '+': 1. / 3., '-': -1. / 3.}
    # get definitive data from KCG's file
    with open(kppdef, 'r') as f:
        ll = f.readlines()
    with open(kppout, 'w') as f:
        for l in ll:
            kparr = l[:32].split()
            kpval = kparr[1:]
            kpdt = datetime.strptime(kparr[0], dformat)
            for t, k in zip(kptbinmid, kpval):
                kp = int(k[0]) + kpdict[k[1]]
                kpmjd = Time(kpdt).mjd + t / 24.
                # kp, kpval, dt, Time(dt,format='mjd').iso
                outstring = f"{kpmjd:.4f} {Time(kpmjd, format='mjd').iso} {k} {kp:.4f}"
                #print(outstring)
                f.write(outstring + "\n")
        # add QL data
        print('Adding QL data')
        req = requests.get(pkpql)
        pkpqldata = req.text.split('\n')
        for l in pkpqldata:
            try:
                kparr = l[:32].split()
                kpval = kparr[1:]
                kpdt = datetime.strptime(kparr[0], dformat)
                for t, k in zip(kptbinmid, kpval):
                    kp = int(k[0]) + kpdict[k[1]]
                    kpmjd = Time(kpdt).mjd + t / 24.
                    # kp, kpval, dt, Time(dt,format='mjd').iso
                    outstring = f"{kpmjd:.4f} {Time(kpmjd, format='mjd').iso} {k} {kp:.4f}"
                    print(outstring)
                    f.write(outstring + "\n")
            except Exception as errmsg:
                print(errmsg)

if __name__ == '__main__':
    outfile = '/software/github/nicergof/utils/bin/tmp/kp.fits'
    caldir = "/software/github/nicergof/utils/bin/tmp"
    dgddir = "/software/github/nicergof/utils/bin/tmp/DGDdir"
    ftpsite = "ftp.swpc.noaa.gov"
    rdir = "/pub/indices/old_indices"
    mk_kp_fits(outfile=outfile, header=None, DGDdir=dgddir,
               getdata=True, verbose=True, clobber=True,
               noaoftp=ftpsite, noaoftpdir=rdir)

