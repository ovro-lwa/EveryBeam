#!/usr/bin/env python3
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

""" Add beam info to OSKAR generated MeasurementSet """


import casacore.tables as pt
import pandas as pd
from lofarantpos import geo
from numpy.linalg import norm
import numpy as np
import logging
import sys


def read_telescope_center(oskar_telescope_dir: str):
    """Read telescope center from OSKAR metadata and convert to ITRF"""
    oskar_telescope_center = pd.read_csv(f"{oskar_telescope_dir}/position.txt",
                                         header=None, sep=" ")

    if len(oskar_telescope_center.columns) == 2:
        oskar_telescope_center["Height"] = [0.]

    oskar_telescope_center.columns = ["Lon", "Lat", "Height"]
    telescope_center_itrf = geo.xyz_from_geographic(np.deg2rad(oskar_telescope_center["Lon"][0]),
                                                    np.deg2rad(oskar_telescope_center["Lat"][0]),
                                                    oskar_telescope_center["Height"][0])
    return telescope_center_itrf


def fix_antenna(oskar_ms_name: str, telescope_center_itrf: np.array):
    """Make POSITION in ::ANTENNA subtable absolute"""
    anttable = pt.table(f"{oskar_ms_name}::ANTENNA", readonly=False, ack=False)
    if np.isclose(norm(anttable[0]["POSITION"]), 0, atol=1e-3):
        anttable.putcol("POSITION", anttable.getcol("POSITION") + telescope_center_itrf)
    else:
        logging.info("Positions are already absolute")
    anttable.close()


def add_array_center(oskar_ms_name: str, telescope_center_itrf: np.array):
    """Add ARRAY_CENTER column to ::OBSERVATION subtable"""
    obstable = pt.table(f"{oskar_ms_name}::OBSERVATION", readonly=False, ack=False)
    if 'ARRAY_CENTER' not in obstable.colnames():
        anttable = pt.table(f"{oskar_ms_name}::ANTENNA", ack=False)
        coldesc = anttable.getcoldesc("POSITION")
        coldesc["name"] = "ARRAY_CENTER"
        coldesc["comment"] = "Reference position for array"
        obstable.addcols(coldesc)

    obstable.putcol("ARRAY_CENTER", np.array([telescope_center_itrf]))
    obstable.close()


def add_phased_array_table(oskar_ms_name: str):
    """Add PHASED_ARRAY subtable to measurement set"""

    anttable = pt.table(f"{oskar_ms_name}::ANTENNA", ack=False)
    phasedarraytable = pt.table(f"{oskar_ms_name}/PHASED_ARRAY",
                                pt.maketabdesc([]),
                                nrow=anttable.nrows())

    oskar_ms = pt.table(f"{oskar_ms_name}", readonly=False, ack=False)
    oskar_ms.putkeyword("PHASED_ARRAY", phasedarraytable, makesubrecord=True)
    oskar_ms.close()

    position_coldesc = anttable.getcoldesc("POSITION")
    position_coldesc['comment'] = 'Position of antenna ﬁeld'
    position_coldesc['name'] = 'POSITION'
    phasedarraytable.addcols(position_coldesc)

    coordinate_system_coldesc = pt.makearrcoldesc("COORDINATE_AXES", 0.,
                                                  shape=[3, 3], comment="Local coordinate system",
                                                  valuetype='double',
                                                  keywords={'QuantumUnits': ['m', 'm', 'm'],
                                                            'MEASINFO': {'Ref': 'ITRF', 'type': 'direction'}})
    phasedarraytable.addcols(coordinate_system_coldesc)
    pt.taql("UPDATE $phasedarraytable SET COORDINATE_AXES=0.");

    element_offset_coldesc = pt.makearrcoldesc("ELEMENT_OFFSET", 0.,
                                               ndim=2, comment="Offset per element",
                                               valuetype='double',
                                               keywords={'QuantumUnits': ['m', 'm', 'm'],
                                                         "MEASINFO": {"type": "position",
                                                                      "Ref": "ITRF"}})
    phasedarraytable.addcols(element_offset_coldesc)

    element_flag_coldesc = pt.makearrcoldesc("ELEMENT_FLAG", 0.,
                                             ndim=2, comment="Offset per element",
                                             valuetype='bool')
    phasedarraytable.addcols(element_flag_coldesc)
    phasedarraytable.close()


def fill_phased_array(oskar_ms_name: str, oskar_telescope_dir: str):
    """Fill the ::PHASED_ARRAY subtable with info from the OSKAR directory"""
    element_locations = pd.read_csv(f"{oskar_telescope_dir}/station/layout.txt",
                                    header=None, sep=" ")
    if len(element_locations.columns) == 2:
        element_locations["Z"] = np.zeros_like(element_locations[0])

    oskar_telescope_center = pd.read_csv(f"{oskar_telescope_dir}/position.txt",
                                         header=None, sep=" ")
    oskar_telescope_center.columns = ["Lon", "Lat"]
    telescope_center_itrf = read_telescope_center(oskar_telescope_dir)

    normal_vector_ellipsoid = geo.normal_vector_ellipsoid(np.deg2rad(oskar_telescope_center["Lon"][0]),
                                                          np.deg2rad(oskar_telescope_center["Lat"][0]))

    local_to_itrf_projection_matrix = geo.projection_matrix(telescope_center_itrf,
                                                            normal_vector_ellipsoid)

    element_locations_itrf = local_to_itrf_projection_matrix @ element_locations.to_numpy().T
    nr_pol = 2
    all_unflagged = np.zeros((nr_pol, element_locations_itrf.shape[1]))

    phasedarraytable = pt.table(f"{oskar_ms_name}::PHASED_ARRAY", readonly=False, ack=False)

    for stationnr in range(len(phasedarraytable)):
        phasedarraytable.putcell("COORDINATE_AXES", stationnr, local_to_itrf_projection_matrix.T)
        phasedarraytable.putcell("ELEMENT_OFFSET", stationnr, element_locations_itrf)
        phasedarraytable.putcell("ELEMENT_FLAG", stationnr, all_unflagged)

    anttable = pt.table(f"{oskar_ms_name}::ANTENNA", ack=False)
    phasedarraytable.putcol("POSITION", anttable.getcol("POSITION"))
    phasedarraytable.close()


def main(oskar_ms_name: str, oskar_telescope_dir: str):
    """Add beam info to OSKAR generated MeasurementSet"""
    telescope_center_itrf = read_telescope_center(oskar_telescope_dir)
    fix_antenna(oskar_ms_name, telescope_center_itrf)
    add_array_center(oskar_ms_name, telescope_center_itrf)
    add_phased_array_table(oskar_ms_name)
    fill_phased_array(oskar_ms_name, oskar_telescope_dir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
