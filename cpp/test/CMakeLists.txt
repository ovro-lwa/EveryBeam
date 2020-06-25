#------------------------------------------------------------------------------
# Add executable(s)
add_executable(tstation tstation.cc)
target_link_libraries(tstation PUBLIC everybeam)
target_link_libraries(tstation PUBLIC OpenMP::OpenMP_CXX)

#------------------------------------------------------------------------------
# Add test
add_test(station-tests tstation)