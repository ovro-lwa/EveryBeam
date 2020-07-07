#------------------------------------------------------------------------------
# CMake file for compiling a comparison between OSKAR and EveryBeam
add_executable(comparison-oskar-generate-beampattern main.cpp)

target_link_libraries(comparison-oskar-generate-beampattern oskar)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/position.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/station000/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/station000/telescope.tm")

#------------------------------------------------------------------------------
# comparison-oskar knits together the cpp code and the python scripts     
add_custom_target(comparison-oskar
     COMMAND ${CMAKE_COMMAND} -E env EXTRA_PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts"
     "${CMAKE_CURRENT_SOURCE_DIR}/generate_basefunction_plots.sh"
     DEPENDS comparison-oskar-generate-beampattern
)

add_test(NAME comparison-oskar-test
     CONFIGURATIONS integration
     COMMAND make comparison-oskar)

