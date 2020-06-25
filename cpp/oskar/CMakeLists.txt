#------------------------------------------------------------------------------
# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------
# build liboskar.so
add_library(oskar SHARED
    OSKARElementResponse.cc
    OSKARDatafile.cc
    OSKARDataset.cc
    oskar_evaluate_dipole_pattern.cc
    oskar_evaluate_spherical_wave_sum.cc)

string(TOLOWER ${CMAKE_PROJECT_NAME} projectname )
set_target_properties(oskar PROPERTIES LIBRARY_OUTPUT_NAME "${projectname}-oskar")

#------------------------------------------------------------------------------
# Link against HDF5 and OpenMP
target_link_libraries(oskar ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})

#------------------------------------------------------------------------------
# install libeverybeam-oskar.so
install(
  TARGETS oskar
  EXPORT EveryBeamTargets
  DESTINATION lib)
