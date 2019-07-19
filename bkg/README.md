This package contains functions which can be used to create an estimated NICER background spectrum based on the
"environmental" model developed by the NICER Guest Observer Facility. The environmental model uses a combination of
the cut-off rigidity (COR_SAX) and the Planetary K index (KP) which gives an estimate of the
space weather environment.  This model also uses the SUN_ANGLE parameter which helps describe
the low-energy background produced by optical loading.  COR_SAX and SUN_ANGLE are contained in the
"make filter" file (either the standard auxil/ni*.mkf file distributed with processed data or
-recommended- the augmented MKF file produced by the "niprefilter2" tool distributed with the NICERDAS HEASoft package).

The KP values are not currently included in either the .mkf or the enhanced (.mkf2) makefilter files,
and must be added to the .mkf2 (or .mkf) file using the "add_kp" function defined here.

Note that this is a PYTHON 3 package

The package also includes a unit test included as an example and a jupyter notebook showing an example of
usage in the testdata subdirectory (which also includes a test .pha file and the auxil
directory from the relevant observation) .


INSTALLATION:
    * Download the bkg_estimator.tar file
    * Untar it in and place the nicergof directory in your $PYTHONPATH
    * Import it  as
            >>> from nicergof.bkg import bkg_estimator as bk
    * There are command line executables, niaddkp and nibkgestimator, which are located in the nicergof/bkg/bin
    directory. You can make an estimated background using these executables from the command line
    by putting them in one of your $PATH directories and making sure that you have execute privilege:
        % chmod u+x nibkgestimator


This background creation method uses two files:
    a) The background events file 30nov18targskc_enhanced.evt
    (current version: https://heasarc.gsfc.nasa.gov/FTP/caldb/data/nicer/xti/pcf/30nov18targskc_enhanced.evt)
    b) the KP.fits file
    (the current version, updated daily, is at https://heasarc.gsfc.nasa.gov/FTP/caldb/data/gen/pcf/kp.fits)

You can access these files virtually by specifying the URLs given above (the default for the functions defined below)
or you can download them to a local directory for faster access.

USAGE EXAMPLES:

Assumptions:
      - the NICER data (a pha file called "test.pha" and an mkf2 file called "1200040103/auxil/ni1200040103.mkf2") is stored in subdirectories of a directory called "testdata" of the current working directory.

    EXAMPLE 1: Command-line execution

        After placing the executables niaddkp and nibkgestimator in your $PATH with execcute permission,
        you can create a nicer background

        a) create the mkf3 file:
            % niaddkp testdata/1200040103/auxil/ni1200040103.mkf2

        b) create the background spectrum:
            % nibkgestimator testdata/test.pha testdata/1200040103/auxil/ni1200040103.mkf3


    EXAMPLE 2: Python usage

        1) import the package:

            >>> from nicergof.bkg import bkg_estimator as be

        2) update the .mkf2 file to include the KP values:

            >>> status = be.add_kp("testdata/1200040103/auxil/ni1200040103.mkf2")

        This will create a ".mkf3" file in the 11200040103/auxil directory

        3) use the mk_bkg_spec_evt function to create the background spectrum:

            >>> bkg_chan, bkgspectot, btotexpo = be.mk_bkg_spec_evt('testdata/test.pha', mkf3file='testdata/1200040103/auxil/ni1200040103.mkf')

        This will create a HEASARC-compliant background PHA file (with .pha replaced by _bkg.pha, i.e.
        the background file for "./test.pha" is "./test_bkg.pha").  It will also return the background channels,
        spectrum and exposure as the python arrays bkg_chan, bkgspectot, btotexpo, which can be plotted or manipulated
        using standard python tools.

CAVEATS:

    * This is PRE-RELEASE software.
    * This software estimates the instrumental background.  X-ray cosmic background (sky background) appropriate to your source is NOT included in the estimated background spectrum produced by mk_bkg_spec_evt().  However, cosmic X-ray background in the NICER blank fields is included in the estimated background.
    * Note that the background events file excludes FPMs #14 & 34, so it uses data from 50 out of the 52 active FPMs.  The estimated background
    * There may be combinations of (KP, COR_SAX, SUN_ANGLE) which are not contained in the background events file; in this case these times are ignored in the output background spectrum.
    * there may be other parameters that are important in determining background, or other parameters which give a better estimate of the background.  This is still under investigation.
    * there are undoubtedly other issues.


Please send any usage issues or problems with the backgrounds produced to nicerhelp@bigbang.gsfc.nasa.gov
