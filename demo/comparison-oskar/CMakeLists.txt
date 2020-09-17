#------------------------------------------------------------------------------
# CMake file for compiling a comparison between OSKAR and EveryBeam
add_executable(make_element_response_image make_element_response_image.cpp)
target_link_libraries(make_element_response_image oskar)

add_executable(make_station_response_image make_station_response_image.cpp)
target_link_libraries(make_station_response_image everybeam)

# Required to get the config.h header
target_include_directories(make_station_response_image PRIVATE "${CMAKE_BINARY_DIR}")


file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/position.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/telescope.tm")

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/telescope.tm/station000/layout.txt"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/station000/telescope.tm")

#------------------------------------------------------------------------------
# comparison-oskar knits together the cpp code and the python scripts     
add_custom_target(comparison-oskar-basefunctions
     COMMAND ${CMAKE_COMMAND} -E env
         EXTRA_PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts"
     "${CMAKE_CURRENT_SOURCE_DIR}/generate_basefunction_plots.sh"
     DEPENDS make_element_response_image
)

add_custom_target(comparison-oskar-station-response
     COMMAND ${CMAKE_COMMAND} -E env
         EXTRA_PATH="${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_SOURCE_DIR}/scripts/coeff_scripts:${CMAKE_SOURCE_DIR}/scripts/misc"
         DATA_DIR=${CMAKE_SOURCE_DIR}/test_data
     "${CMAKE_CURRENT_SOURCE_DIR}/compare_stationresponse.sh"
     DEPENDS make_station_response_image
)

