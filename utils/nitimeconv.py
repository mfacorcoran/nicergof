#!/usr/bin/env python

__name__    = 'nitimeconv'
__author__  = 'Teruaki Enoto'
__version__ = '1.01'
__date__    = '2018 April 7'

from optparse import OptionParser
from astropy.time import Time 

NICER_MJDREFI   = 56658.0
NICER_MJDREFF   = 0.000777592592592593
NICER_TIMEZERO  = 0.0
LEAP_INIT = 2.0

NICER_MET_ORIGIN   = Time('2014-01-01T00:00:00.000',format='isot',scale='utc')
FERMI_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
NUSTAR_MET_ORIGIN  = Time('2010-01-01T00:00:00.000',format='isot',scale='utc')
RXTE_MET_ORIGIN    = Time('1994-01-01T00:00:00.000',format='isot',scale='utc')
SUZAKU_MET_ORIGIN  = Time('2000-01-01T00:00:00.000',format='isot',scale='utc')
SWIFT_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
XMM_MET_ORIGIN     = Time('1998-01-01T00:00:00.000',format='isot',scale='tt')
CHANDRA_MET_ORIGIN = XMM_MET_ORIGIN
"""
Fermi seconds since 2001.0 UTC (decimal)		410227203.000
Fermi mission week (integer)		291
LIGO/GPS seconds since 1980-01-06 UTC (decimal)		1072569616.000
NuSTAR seconds since 2010.0 UTC (decimal)		126230401.000
RXTE seconds since 1994.0 UTC (decimal)		631152007.000
RXTE seconds since 1994.0 UTC (hexadecimal)		0x259e9d87
RXTE mission day number (dddd:hh:mm:ss)		7305:00:00:07.000
RXTE decimal mission day (dddd.ddd...)		7305.00008102
Suzaku seconds since 2000.0 UTC (decimal)		441849603.000
Swift seconds since 2001.0 UTC (decimal)		410227211.463
XMM/Chandra seconds since 1998.0 TT (decimal)		504921667.184
"""

usage = """ 
NAME 
	nitimeconv - Convert NICER related time into different time systems.

USAGE 
	%prog intime -f format -s scale 

DESCRIPTION
	'%prog' takes an input time with specified format and scale, and then
	converts it into another time systems. This gives similar outputs of 
	the HEASARC "xTime - A Date/Time Conversion Utility" 

		https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl

	The NICER Mission Elapse Time (MET), the "TIME" column is defined as 
	elapsed TT seconds since the epoch 2014-01-01T00:00:00 UTC. Conversion 
	of this NICER timestamps to absolute time in TT_MJD can be defined as 

		MJD(TT) = (MJDREFI+MJDREFF) + (TIMEZERO+TIME)/86400
		MJD(UTC) = (MJDREFI) + (TIMEZERO+TIME+LEAPINIT=2)/86400

	This script converts the MET to MJD(TT) using the above definition, and
	then calculates MJD(UTC) via the python library 'astropy' witout using 
	the above second equation. It is recommended to avoid representing UTC 
	in MJD unit. The NICER MET is compared with other METs of different 
	X-ray missions (Chandra/XMM, RXTE, Fermi, and Suzaku). The MET sometime
	requires further time corrections, and please double check the xTime. 

REFERENCES
	https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html

	http://docs.astropy.org/en/stable/api/astropy.time.Time.html#astropy.time.Time

EXAMPLES
	1. Get the NICER TIME origion. 
		$ %prog 2014-01-01T00:00:00 -f isot -s utc

	2. Convert Calender format to the NICER time (126266402.000).
		$ %prog 2018-01-01T10:00:00 -f isot -s utc

	3. Convert NICER TIME to other time format (oposite of 2).
		$ %prog 126266402 -f met -s met 

	4. Convert MJD_UTC 58000 to other formats. 
		$ %prog 58000 -f mjd -s utc  

	5. Convert Fermi MET to other formats.
		$ %prog 13000 -f met -s met_fermi

	6. Convert elapsed day format (UTC) to other formats. 
		nitimeconv.py 2017:051:03:20:00.000 -f yday -s utc

"""
parser = OptionParser(usage=usage)
parser.add_option("-f","--format",dest="format",default="met",
       action="store",help="Time format, any of (met, jd, mjd, isot, and yday)",type="string")
parser.add_option("-s","--scale",dest="scale",default="met",
       action="store",help="Time scale, any of (utc, tt, met, met_nicer, met_fermi, met_nustar, met_suzaku, met_xmm, and met_chandra)",type="string")
(options, args) = parser.parse_args()

if len(args) != 1:
	print("try: %s.py --help" % __name__)
	print("usage: %s.py intime -f fomrat -s scale" % __name__)	
	quit()
input_value = args[0]

dump  = "----- Input Time Value and Formats -----\n"
dump += "intime: %s\n" % str(input_value)
dump += "format: %s\n" % options.format
dump += "scale : %s\n" % options.scale 

if options.format == "met":
	if options.scale == "met" or options.scale == "met_nicer":
		mission = float(input_value)
		mjd_tt  = NICER_MJDREFI+NICER_MJDREFF+(NICER_TIMEZERO+mission)/86400.0
		time_tt = Time(mjd_tt,format='mjd',scale='tt')
		time_utc = time_tt.utc
	elif options.scale == "met_fermi":
		time_utc = Time(float(input_value)+FERMI_MET_ORIGIN.gps,format='gps',scale='utc')
		time_tt  = time_utc.tt
		mission  = time_tt.gps + NICER_MET_ORIGIN.gps
	elif options.scale ==  "met_nustar":
		time_utc = Time(float(input_value)+NUSTAR_MET_ORIGIN.gps,format='gps',scale='utc')
		time_tt  = time_utc.tt
		mission  = time_tt.gps + NICER_MET_ORIGIN.gps	
	elif options.scale == "met_suzaku":
		time_utc = Time(float(input_value)+SUZAKU_MET_ORIGIN.gps,format='gps',scale='utc')
		time_tt  = time_utc.tt
		mission  = time_tt.gps + NICER_MET_ORIGIN.gps		
	elif options.scale == "met_xmm" or options.scale == "met_chandra":
		time_utc = Time(float(input_value)+XMM_MET_ORIGIN.gps,format='gps',scale='utc')
		time_tt  = time_utc.tt
		mission  = time_tt.gps + NICER_MET_ORIGIN.gps	
else:
	if input_value.isdigit():
		time = Time(float(input_value),format=options.format,scale=options.scale)
	else:
		time = Time(str(input_value),format=options.format,scale=options.scale)		
	time_tt = time.tt 
	time_tt.format = 'mjd'
	mjd_tt = time_tt 
	mission = (float(mjd_tt.mjd) - NICER_MJDREFI - NICER_MJDREFF) * 86400.0 - NICER_TIMEZERO
	time_utc = time.utc 
	time_utc.format = 'mjd'

fermi_time  = time_tt.gps - FERMI_MET_ORIGIN.gps
nustar_time = time_tt.gps - NUSTAR_MET_ORIGIN.gps 
#rxte_time   = time_tt.gps - RXTE_MET_ORIGIN.gps 
suzaku_time = time_tt.gps - SUZAKU_MET_ORIGIN.gps 
#swift_time  = time_tt.gps - SWIFT_MET_ORIGIN.gps 
xmm_time    = time_tt.gps - XMM_MET_ORIGIN.gps 
chandra_time = time_tt.cxcsec

dump += "----- Calendar Time Formats -----\n"
dump += "ISO8601_TT : %s (TT)\n" % time_tt.isot 
dump += " JD_TT     : %.8f (TT) \n" % time_tt.jd
dump += "MJD_TT     :   %.8f (TT)\n" % time_tt.mjd
dump += "ISO8601_UTC: %s (UTC)\n" % time_utc.isot 
dump += " JD_UTC    : %.8f (UTC) \n" % time_utc.jd
dump += "MJD_UTC    :   %.8f (UTC) \n" % time_utc.mjd
dump += "----- Mission-Specific Time Formats (Misson Elapsed Time, NET) -----\n"
dump += "Fermi seconds sicne 2001.0 UTC (decimal)     : %.6f\n" % fermi_time
dump += "NuSTAR seconds since 2010.0 UTC (decimal)    : %.6f\n" % nustar_time
#dump += "RXTE seconds since 1994.0 UTC (decimal)      : %.8f\n" % rxte_time 
dump += "Suzaku seconds since 2000.0 UTC (decimal)    : %.6f\n" % suzaku_time
#dump += "Swift seconds since 2001.0 UTC (decimal): %.8f\n" % swift_time
dump += "XMM seconds since 1998.0 TT (decimal)        : %.6f\n" % xmm_time
dump += "Chandra seconds since 1998.0 TT (decimal)    : %.6f\n" % chandra_time
dump += "NICER seconds since 2014.0 UTC (decimal)     : %.6f" % mission 
print(dump)



