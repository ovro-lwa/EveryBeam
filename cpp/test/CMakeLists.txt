#------------------------------------------------------------------------------
# Add executable(s)
add_executable(tstation tstation.cc)
target_link_libraries(tstation PUBLIC everybeam)
target_link_libraries(tstation PUBLIC OpenMP::OpenMP_CXX)

add_executable(taocommon taocommon.cc)
add_executable(teigen teigen.cc)
add_executable(tpybind11 tpybind11.cc)
target_link_libraries(tpybind11 PRIVATE pybind11::embed)
#------------------------------------------------------------------------------
# Add test
add_test(station-tests tstation)
add_test(submodule-tests taocommon)
add_test(submodule-tests teigen)
add_test(submodule-tests tpybind11)