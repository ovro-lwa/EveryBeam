find_package(Doxygen)

if(NOT DOXYGEN_FOUND)
  add_custom_target(
    doxygen
    COMMAND
      echo -e
      "\\nDoxygen is not found. Please install Doxygen, rerun cmake and try again.\\n"
    COMMAND exit 1
    VERBATIM)
else()
  set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

  #Replace variables inside @@ with the current values
  configure_file(Doxyfile.in Doxyfile @ONLY)

  add_custom_target(
    doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ../Doxyfile
    WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR}
    COMMENT "Generating documentation with Doxygen")
endif()

find_package(Sphinx)

if(NOT SPHINX_FOUND)
  add_custom_target(
    doc
    COMMAND
      echo -e
      "\\nSphinx is not found. Please install Sphinx, rerun cmake and try again.\\n"
    COMMAND exit 1
    VERBATIM)
else()
  set(EVERYBEAM_SO_PATH ${CMAKE_BINARY_DIR}/python)
  message(${EVERYBEAM_SO_PATH})

  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/html)

  add_custom_target(
    doc
    ${CMAKE_COMMAND}
    -E
    env
    EVERYBEAM_SO_PATH=${EVERYBEAM_SO_PATH}
    ${SPHINX_EXECUTABLE}
    -b
    html
    ${SPHINX_SOURCE}
    ${SPHINX_BUILD}
    # Tell Breathe where to find the Doxygen output
    -Dbreathe_projects.EveryBeam=${DOXYGEN_OUTPUT_DIR}/xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx"
    DEPENDS doxygen)
endif()
