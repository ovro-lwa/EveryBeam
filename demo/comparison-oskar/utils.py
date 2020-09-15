#!/usr/bin/env python3
"""
Common functions used in the demo project.
"""

import sys
import numpy as np
from astropy.time import Time, TimeDelta
import oskar

"""
Checks if the maximum difference between two matrices are no bigger than the tolerance.
"""
def check_tolerance(tolerance, A, B):
    if tolerance:
        difference = np.nanmax(np.abs(A - B))
        if difference > tolerance:
            sys.exit(
                "Difference between OSKAR and EveryBeam spherical wave model is {}, which is larger than the tolerance {}".format(
                    difference, tolerance
                )
            )


"""
Create a settings file.
"""
def create_settings(app, settings_path, current_settings):
    open(settings_path, 'w').close()
    settings = oskar.SettingsTree(app, settings_path)
    settings.from_dict(current_settings)
    return settings


"""
Returns optimal start time for field RA and observation length.
"""
def get_start_time(ra0_deg, length_sec, optimal_time):
    t = Time(optimal_time, scale='utc', location=('116.764d', '0d'))
    dt_hours = (24.0 - t.sidereal_time('apparent').hour) / 1.0027379
    dt_hours += (ra0_deg / 15.0)
    start = t + TimeDelta(dt_hours * 3600.0 - length_sec / 2.0, format='sec')
    return start.value