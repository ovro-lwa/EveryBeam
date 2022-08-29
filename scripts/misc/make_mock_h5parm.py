#!/usr/bin/env python3
from idg.h5parmwriter import H5ParmWriter
import numpy as np
from datetime import datetime

"""
Script for generating the MOCK_H5PARM.h5 file, used for unittesting the
H5ParmATerm implementation.
Please note that running this script requires the python idg module
to be on your PYTHONPATH.
"""

## INPUT
fname = "MOCK_H5PARM.h5"
solsetname = "coefficients000"

antenna_names = np.array(["Antenna0", "Antenna1"])
antenna_positions = np.array([[0.0, 1.0, 2.0], [3, 4, 5]])

nr_stations = antenna_names.size
# Amplitude --> Second order polynomial
nr_parameters_ampl = 6
time_array_ampl = np.linspace(0, 10, num=3)

# Phase --> First order polynomial
nr_parameters_phase = 3
time_array_phase = np.linspace(0, 10, num=11)

image_size = 1.0
subgrid_size = 32
## END INPUT


def init_h5parm_solution_table(
    h5parm_object,
    soltab_type,
    axes_info,
    antenna_names,
    time_array,
    image_size,
    subgrid_size,
    basisfunction_type="lagrange",
    history="",
):
    soltab_info = {
        "amplitude": "amplitude_coefficients",
        "phase": "phase_coefficients",
    }

    assert soltab_type in soltab_info.keys()
    soltab_name = soltab_info[soltab_type]

    h5parm_object.create_solution_table(
        soltab_name,
        soltab_type,
        axes_info,
        dtype=np.float_,
        history=f'CREATED at {datetime.today().strftime("%Y/%m/%d")}; {history}',
    )

    # Set info for the "ant" axis
    h5parm_object.create_axis_meta_data(
        soltab_name, "ant", meta_data=antenna_names
    )

    # Set info for the "dir" axis
    h5parm_object.create_axis_meta_data(
        soltab_name,
        "dir",
        attributes={
            "basisfunction_type": basisfunction_type,
            "image_size": image_size,
            "subgrid_size": subgrid_size,
        },
    )

    # Set info for the "time" axis
    h5parm_object.create_axis_meta_data(
        soltab_name, "time", meta_data=time_array
    )
    return h5parm_object


axes_labels = ["ant", "time", "dir"]
axes_data_amplitude = dict(
    zip(
        axes_labels,
        (nr_stations, time_array_ampl.size, nr_parameters_ampl),
    )
)

axes_data_phase = dict(
    zip(
        axes_labels,
        (
            nr_stations,
            time_array_phase.size,
            nr_parameters_phase,
        ),
    )
)

h5writer = H5ParmWriter(fname, solution_set_name=solsetname, overwrite=True)


h5writer.add_antennas(antenna_names, antenna_positions)

h5writer = init_h5parm_solution_table(
    h5writer,
    "amplitude",
    axes_data_amplitude,
    antenna_names,
    time_array_ampl,
    image_size,
    subgrid_size,
)

h5writer = init_h5parm_solution_table(
    h5writer,
    "phase",
    axes_data_phase,
    antenna_names,
    time_array_phase,
    image_size,
    subgrid_size,
)


amplitude_coeffs = np.array([[1.0, 2.0, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]])
for i, tampl in enumerate(time_array_ampl):
    amplitude_coeffs_tmp = amplitude_coeffs * (i + 1)
    amplitude_coeffs_tmp = amplitude_coeffs_tmp.reshape(
        nr_stations, 1, nr_parameters_ampl
    )
    offset_amplitude = (0, i, 0)
    h5writer.fill_solution_table(
        "amplitude_coefficients", amplitude_coeffs_tmp, offset_amplitude
    )

phase_coeffs = np.array([[1, 2, 3], [4, 5, 6]])
for i, tphase in enumerate(time_array_phase):
    phase_coeffs_tmp = phase_coeffs * (i + 1)
    phase_coeffs_tmp = phase_coeffs_tmp.reshape(
        nr_stations, 1, nr_parameters_phase
    )
    offset_phase = (0, i, 0)
    h5writer.fill_solution_table(
        "phase_coefficients", phase_coeffs_tmp, offset_phase
    )
