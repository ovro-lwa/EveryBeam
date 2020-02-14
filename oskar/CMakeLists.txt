# directory for config.h
include_directories(${CMAKE_BINARY_DIR})

# build liboskar.so
add_library(oskar SHARED
    OSKARElementResponse.cc
    OSKARSphericalWaveCoeff.cc
    oskar_evaluate_dipole_pattern.cc
    oskar_evaluate_spherical_wave_sum.cc)

target_link_libraries(oskar ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})

# install liboskar.so
install(TARGETS oskar DESTINATION lib)
