# Find boost
find_package(Boost COMPONENTS unit_test_framework)

add_executable(taocommon taocommon.cc)
add_executable(teigen teigen.cc)
add_executable(tpybind11 tpybind11.cc)
target_link_libraries(tpybind11 PRIVATE pybind11::embed)

# Add test
add_test(submodule-tests taocommon)
add_test(submodule-tests teigen)
add_test(submodule-tests tpybind11)

#------------------------------------------------------------------------------
# Add tests if Boost found
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    set(TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/runtests.cc)

    # Test LOFAR beamforming
    if( LOFAR_MOCK_MS )        
        set(TEST_SOURCES 
            ${TEST_SOURCES}
            ${CMAKE_CURRENT_SOURCE_DIR}/tload_lofar.cc
        )
    endif()

    # Test VLA beamforming
    if( VLA_MOCK_MS )
        set(TEST_SOURCES 
            ${TEST_SOURCES}
            ${CMAKE_CURRENT_SOURCE_DIR}/tvla.cc
        )
    endif()

    # Test MWA beamforming (assumes that mwa h5 coefficients file can be found!)
    if( MWA_MOCK_MS AND MWA_COEFF_PATH )
        set(TEST_SOURCES 
            ${TEST_SOURCES}
            ${CMAKE_CURRENT_SOURCE_DIR}/tmwa.cc
        )
    elseif( MWA_MOCK_MS OR MWA_COEFF_PATH )
        message(WARNING "To compile the MWA test, both the MWA_MOCK_MS and MWA_COEFF_PATH needs to be set. Now skipping MWA test compiling.")
    endif()

    # Add test executable
    add_executable(runtests ${TEST_SOURCES})
    target_link_libraries(runtests everybeam ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    # Required to get the config.h header
    target_include_directories(runtests PRIVATE "${CMAKE_BINARY_DIR}")

    add_test(runtests runtests)
endif()


