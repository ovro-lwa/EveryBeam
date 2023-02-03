from astropy.coordinates import AltAz, EarthLocation, ITRS, SkyCoord
from astropy.time import Time

import astropy.units as u
import casacore.tables as ct
import everybeam as eb
import numpy as np


def radec_to_xyz(ra, dec, time):
    """Convert ra and dec coordinates to ITRS coordinates for LOFAR observations.

    Args:
        ra (astropy Quantity): right ascension
        dec (astropy Quantity): declination
        time (float): MJD time in seconds
    Returns:
        pointing_xyz (ndarray): NumPy array containing the X, Y and Z coordinates
    """
    obstime = Time(time / 3600 / 24, scale="utc", format="mjd")
    loc_LOFAR = EarthLocation(
        lon=0.11990128407256424,
        lat=0.9203091252660295,
        height=6364618.852935438 * u.m,
    )

    dir_pointing = SkyCoord(ra, dec)
    dir_pointing_altaz = dir_pointing.transform_to(
        AltAz(obstime=obstime, location=loc_LOFAR)
    )
    dir_pointing_xyz = dir_pointing_altaz.transform_to(ITRS)

    pointing_xyz = np.asarray(
        [dir_pointing_xyz.x, dir_pointing_xyz.y, dir_pointing_xyz.z]
    )
    return pointing_xyz


# Set path to LOFAR LBA MS and load telescope
ms_path = "/path/to/my.ms"

telescope = eb.load_telescope(ms_path)
assert type(telescope) == eb.LOFAR

# Time slots at which to evaluate the beam response.
ms_times = ct.taql("SELECT TIME FROM {ms:s}".format(ms=ms_path))
times = ms_times.getcol("TIME")[0:5]

# Frequencies at which to evaluate the beam response.
ms_freqs = ct.taql(
    "SELECT CHAN_FREQ FROM {ms:s}::SPECTRAL_WINDOW".format(ms=ms_path)
)
freqs = ms_freqs.getcol("CHAN_FREQ").squeeze()

# Obtain the reference direction from the Measurement Set.
ms_dirs = ct.taql("SELECT REFERENCE_DIR FROM {ms:s}::FIELD".format(ms=ms_path))
ra_ref, dec_ref = ms_dirs.getcol("REFERENCE_DIR").squeeze()

# Change these to the RA and DEC of interest, in units of radians.
ra, dec = 1.34, 1.56

# Here we obtain an array of (x, y, z) coordinates for each time slot in the Measurement Set.
# ITRF coordinates of the reference direction.
reference_xyz = list(
    zip(*radec_to_xyz(ra_ref * u.rad, dec_ref * u.rad, times))
)
# ITRF coordinates of the phase centre to correct the array factor for.
phase_xyz = list(zip(*radec_to_xyz(ra * u.rad, dec * u.rad, times)))
# Get teh array factor for station 0
station_id = 0
array_factor = telescope.array_factor(
    times[0], station_id, freqs[0], phase_xyz[0], reference_xyz[0]
)

print(array_factor)
