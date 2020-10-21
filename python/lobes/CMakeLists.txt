project(lobes_pybind11)

# Create the binding library
pybind11_add_module(pylobes 
                    lobes.cc
                    )


target_include_directories(pylobes PUBLIC ${HDF5_LIBRARIES})
target_link_libraries(pylobes PUBLIC ${HDF5_LIBRARIES})
