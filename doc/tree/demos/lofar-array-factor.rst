.. _lofardemoarrayfactor:

LOFAR telescope: correcting for the array factor in a given direction
=====================================================================

The EveryBeam Python bindings can be used to evaluate the response in an arbitrary equatorial direction.
This demo illustrates in a step-by-step manner a minimal example of doing this.
The complete code is shown at the bottom of this page.

It all starts with importing the ``everybeam`` python module, along with any other
libraries that are needed for your project and specifying a path to the Measurement Set
that will be used:

::

    from astropy.coordinates import AltAz, EarthLocation, ITRS, SkyCoord
    from astropy.time import Time

    import astropy.units as u
    import everybeam as eb
    import numpy as np


    # Set path to LOFAR LBA MS and load telescope
    ms_path = "/path/to/my.ms"

The telescope can now be loaded with ``eb.load_telescope``. This should return an instance of a  ``LOFAR`` telescope.

::

    telescope = eb.load_telescope(ms_path)
    assert type(telescope) == eb.LOFAR

Please note that since no optional arguments were passed to ``load_telescope``, no differential beam will applied and the default element response model (``"Hamaker"``) will be used.

In order to evaluate the array_factor, we need to pass some user input. For the array factor, this includes a time, station index, frequency, the direction we want to know the array factor in and the reference direction (where the array factor is unity). For the ``array_factor`` method these coordinates need to be in ITRF coordinates. We therefore also need to convert from equatorial coordinates into these ITRF coordinates.

::

    import casacore.tables as ct
    # Time slots at which to evaluate the beam response.
    ms_times = ct.taql('SELECT TIME FROM {ms:s}'.format(ms=ms_path))
    times = ms_times.getcol('TIME')[0]

    # Frequencies at which to evaluate the beam response.
    ms_freqs = ct.taql('SELECT CHAN_FREQ FROM {ms:s}::SPECTRAL_WINDOW'.format(ms=ms_path))
    freqs = ms_freqs.getcol('CHAN_FREQ').squeeze()

    # Obtain the reference direction from the Measurement Set.
    ms_dirs = ct.taql('SELECT REFERENCE_DIR FROM {ms:s}::FIELD'.format(ms=ms_path))
    ra_ref, dec_ref = dirs.getcol('REFERENCE_DIR').squeeze()

    # Change these to the RA and DEC of interest, in units of radians.
    ra, dec = 1.34, 1.56

Now that we have our initial information we can calculate the ITRF coordinates. We use AstroPy's coordinate transformations to first obtain the local azimuth and elevation of the object as seen by LOFAR's centre. These are then transformed to the ITRS system.

::

    def radec_to_xyz(ra, dec, time):
        ''' Convert ra and dec coordinates to ITRS coordinates for LOFAR observations.
        
        Args:
            ra (astropy Quantity): right ascension
            dec (astropy Quantity): declination
            time (float): MJD time in seconds
        Returns:
            pointing_xyz (ndarray): NumPy array containing the X, Y and Z coordinates
        '''
        obstime = Time(time/3600/24, scale='utc', format='mjd')
        loc_LOFAR = EarthLocation(lon=0.11990128407256424, lat=0.9203091252660295, height=6364618.852935438*u.m)

        dir_pointing = SkyCoord(ra, dec) 
        dir_pointing_altaz = dir_pointing.transform_to(AltAz(obstime=obstime, location=loc_LOFAR))
        dir_pointing_xyz = dir_pointing_altaz.transform_to(ITRS)

        pointing_xyz = np.asarray([dir_pointing_xyz.x, dir_pointing_xyz.y, dir_pointing_xyz.z])
        return pointing_xyz

    # Here we obtain an array of (x, y, z) coordinates for each time slot in the Measurement Set.
    # ITRF coordinates of the reference direction.
    reference_xyz = list(zip(*radec_to_xyz(ra_ref * u.rad, dec_ref * u.rad, times)))
    # ITRF coordinates of the phase centre to correct the array factor for.
    phase_xyz = list(zip(*radec_to_xyz(ra * u.rad, dec * u.rad, times)))


We're now all set to compute the array factor response using the ``array_factor`` method. It does not take arrays as input, so we have to calculate it per time and frequency slot. The array factor for the first time slot would then be:

::

    array_factor = telescope.array_factor(times[0], station_id, freqs[0], phase_xyz[0], reference_xyz[0])

This returns a 2x2 array with the response of the XX, XY, YX and YY. It is recommendable to parallelise this calculation over channels, for example, to speed things along. Since EveryBeam's methods do not take arrays, the full time range can be calculated using something like the following:

::

    # Station to calculate the response for.
    station_id = 23
    # Frequency channel to process.
    ifreq = 0

    # Obtain the array_factor response for the specified station, for the first time and frequency slots.
    freq = freqs[ifreq]
    # The reponse matrix has complex numbers.
    timeslices = np.empty((len(times), 2, 2), dtype=np.complex128)
    for itime, time in enumerate(times):
        if not useElementResponse and useArrayFactor:
            # Array-factor-only correction.
            beam = obs.array_factor(times[itime], stationnum, freq, phase_xyz[itime], reference_xyz[itime])
        else:
            beam = obs.station_response(time=time, station_idx=stationnum, freq=freq, ra=ra, dec=dec)
        timeslices[itime] = beam

A complete overview of the code is shown below:

.. literalinclude:: ../../../demo/python/lofararrayfactor.py
