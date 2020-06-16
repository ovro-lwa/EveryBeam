# This module tries to find the BeamModel library on your system
#
# Once done this will define
#  BEAMMODEL_FOUND       - system has BeamModel
#  BEAMMODEL_INCLUDE_DIR - the BeamModel include directory
#  BEAMMODEL_LIBRARIES   - link these to use the beam model library

find_package(PackageHandleStandardArgs)

find_path(
    BEAMMODEL_INCLUDE_DIR
    NAMES EveryBeam/Station.h
    HINTS ${BEAMMODEL_ROOT_DIR}
    PATH_SUFFIXES include
)

find_library(
    BEAMMODEL_STATION_RESPONSE_LIBRARY
    NAMES everybeam
    HINTS ${BEAMMODEL_ROOT_DIR}
    PATH_SUFFIXES lib
)

set(BEAMMODEL_LIBRARIES CACHE INTERNAL "")
list(APPEND BEAMMODEL_LIBRARIES ${BEAMMODEL_STATION_RESPONSE_LIBRARY})

mark_as_advanced(BEAMMODEL_STATION_RESPONSE_LIBRARY)

find_package_handle_standard_args(BeamModel DEFAULT_MSG BEAMMODEL_LIBRARIES BEAMMODEL_INCLUDE_DIR)
