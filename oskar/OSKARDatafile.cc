#include "OSKARDatafile.h"

#include <iostream>

Datafile::Datafile(
    const std::string& filename)
{
    // Open file
    std::cout << "read oskar datafile: " << filename << std::endl;
    m_h5_file.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

    // Disable HDF5 error prints
    H5::Exception::dontPrint();
};

std::shared_ptr<Dataset> Datafile::get(
    const unsigned int freq)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    // Find dataset for frequency
    auto entry = m_map.find(freq);

    // If found, retrieve pointer to dataset
    if (entry != m_map.end()) {
        return entry->second;
    }

    // Read and return dataset
    std::shared_ptr<Dataset> dataset_ptr;
    dataset_ptr.reset(new Dataset(*m_h5_file, freq));
    m_map.insert({freq, dataset_ptr});
    return dataset_ptr;
}
