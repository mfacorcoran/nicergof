# Using nicer.py

## Introduction

The Neutron star Interior Composition ExploreR, or NICER, is an X-ray observatory installed on the International Space Station in June 2017.  The primary science of NICER is to accurately measure small distortions in the pulse profiles of rapidly spinning neutron stars (X-ray pulsars) caused by the warping of spacetime caused by the neutron stars extreme density, and so constrain the density of the neutron star. NICER also measures and monitors X-ray emission from other cosmic X-ray sources, like black hole systems, AGN, single and binary stars, and other sources.  NICER has the capability of providing high time-resolution studies of changes in the X-ray spectrum of these sources.

NICER data are stored at the [High Energy Astrophysics Science Archive Research Center](https://heasarc.gsfc.nasa.gov) (HEASARC).  More information can be found at the [NICER mission home page](https://heasarc.gsfc.nasa.gov/docs/nicer) at the HEASARC.

NICER data processing and analysis is mainly done using the `nicer` software sub-package distributed with the [`HEASoft`](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/) data analysis software package.

## The nicer.py Python Package

`nicer.py` is a python package to help access and analyze processed nicer data.  Data access and analysis is done by creating a `nicerObs` observation object.

`nicer.py` assumes that the NICER data is located in a directory tree which is the same as that used in the HEASARC archive

`nicerObs` has the following attributes:

* obsid: The NICER observation identification number
* datadir: The directory path to the NICER observation data directory
* log: The directory path to the NICER observation data directory log directory (`datadir/<obsid>/log`)
* aux: The directory path to the NICER observation data directory auxil directory (`datadir/<obsid>/auxil`)
* xti: The directory path to the NICER observation data directory xti data directory (`datadir/<obsid>/xti`)

and the following methods:

* get_eventfile(): return the name of the events FITS file
* get_eventsdf(): get the events list as a pandas DataFrame
* get_prod(): returns a list of data products (in development)
* get_hk(): returns a list of housekeeping files (in development)
* get_gti(): returns the list of gtis (and durations) as a pandas DataFrame
* get_event_flags(): gets the event flag bit values as a numpy array of integers. event flag of 24 = 2**3 + 2**4 (a 1 in bit 3 and 4) corresponds to events detected by both the fast and slow chains (i.e. the "best" events)
* get_event_times(): get the event times from the specified event file
* get_lc():  return the binned lightcurve (bincounts, bincenters and binwidths) for all gtis (or, optionally, selected gtis) in the events file.
* fold(): return a phase-folded lightcurve (in development)
* get_spectrum():  calculate the binned spectrum optionally between tstart and tstart

`nicer.py` also includes the following convenience functions:

* filter_flag(): filter an event list based on the EVENT_FLAGS
* parse_nicer_flag(): translates a NICER bit event flag & prints result
* plot_gti_lc(): plots the lightcurve for the specified obsid for all the gtis
* set_model(): Creates an XSPEC model instance for a two component absorbed thermal model plus gaussian line
* get_obsDF(): Defines a pandas DataFrame summarizing the analysis of a NICER spectrum with blank spectrum parameters
* get_obsid_specparams(): Get spectrum parameters from fit for a NICER spectrum; returns a NICER spectrum DataFrame
* concat_nobsDF(): combines NICER spectrum DataFrames into a single DataFrame
* set_xspec(): sets the pha and model instances for a given obsid; default model is a 2-temperature absorbed vapec model
* fit_spectrum(): fits a spectrum in PyXSPEC using input pha and model instances

# Usage:
Assume that your NICER observation directory is in `$HOME/nicerdata`, and there's an observation directory corresponding to OBSID = 1012010101. The ``nicer`` package should also be in your `$PYTHONPATH`

To create a nicerObs obsid object, you would:

```
ipython> from nicer.nicer import nicerOBS
ipython> import os
ipython> home = os.environ['HOME']
ipython> datadir = os.path.join(home,'nicerdata')
ipython> obsid = 1012010101
ipython> nobs = nicerObs(obsid, rootdir = datadir)
ipython> obsdir = nobs.datadir.replace(home,'$HOME')
ipython> xtidir = nobs.xti.replace(home,'$HOME')
ipython> evtfile = nobs.get_eventfile().replace(home,'$HOME')
ipython> print("XTI data for NICER obsid {obs} is {xti}".format(obs=obsdir, xti=xtidir))
ipython> print("Cleaned event file is {evtf}".format(evtf = evtfile))
ipython> gtis = nobs.get_gti()
ipython> print(gtis)
```
