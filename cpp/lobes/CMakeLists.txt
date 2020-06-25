#------------------------------------------------------------------------------
# Required packages for lobes
# pybind and eigen not required here, will be moved to python dir with pybindings
# find_package(pybind11 REQUIRED)
# find_package (Eigen3 REQUIRED NO_MODULE)

# pybind11_add_module(pylobes SHARED lobes.cc ElementResponse.cc LOBESElementResponse.cc)
# target_link_libraries (pylobes PUBLIC Eigen3::Eigen ${HDF5_LIBRARIES})

add_library(lobes SHARED LOBESElementResponse.cc)
string(TOLOWER ${CMAKE_PROJECT_NAME} projectname )
set_target_properties(lobes PROPERTIES LIBRARY_OUTPUT_NAME "${projectname}-lobes")

#------------------------------------------------------------------------------
# Install
install(
  TARGETS lobes
  EXPORT EveryBeamTargets
  DESTINATION lib)

#------------------------------------------------------------------------------
# TODO: can't we remove all this?
set(MISC_DIR ${CMAKE_SOURCE_DIR}/scripts/misc)
configure_file(${MISC_DIR}/test_lobes.py ${CMAKE_CURRENT_BINARY_DIR}/test_lobes.py COPYONLY)

configure_file(${MISC_DIR}/F4far_new.py ${CMAKE_CURRENT_BINARY_DIR}/F4far_new.py COPYONLY)
configure_file(${MISC_DIR}/calc_modes_1deg_noback.py ${CMAKE_CURRENT_BINARY_DIR}/calc_modes_1deg_noback.py COPYONLY)
configure_file(${MISC_DIR}/hamaker_vs_lobes.py ${CMAKE_CURRENT_BINARY_DIR}/hamaker_vs_lobes.py COPYONLY)
configure_file(${MISC_DIR}/test_beam_model.py ${CMAKE_CURRENT_BINARY_DIR}/test_beam_model.py COPYONLY)

configure_file(${CMAKE_SOURCE_DIR}/coeffs/CS302_coords.mat ${CMAKE_CURRENT_BINARY_DIR}/CS302_coords.mat COPYONLY)

# TODO Too large for git repo, get file from somewhere else
# configure_file(LBA_CS302_fine.mat LBA_CS302_fine.mat COPYONLY)
