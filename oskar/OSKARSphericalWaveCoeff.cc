#include "OSKARSphericalWaveCoeff.h"

Dataset::Dataset(
    H5::H5File& h5_file,
    const unsigned int freq)
{
    // Try to read coefficients for this frequency
    std::string dataset_name = std::to_string((int) (freq / 1e6));

    try {
        // Open dataset
        H5::DataSet dataset = h5_file.openDataSet(dataset_name);

        // Read dataset dimensions
        H5::DataSpace dataspace = dataset.getSpace();
        unsigned int rank = dataspace.getSimpleExtentNdims();
        assert(rank == m_dataset_rank);

        // Get dimensions
        hsize_t dims[rank];
        dataspace.getSimpleExtentDims(dims, NULL);
        m_nr_elements = dims[0];
        assert(dims[1] == m_nr_pols); // pola, polb
        assert(dims[2] == m_nr_tetm); // te, tm
        m_l_max = dims[3]; // l_max
        assert(dims[4] == m_l_max); // m_abs

        // Read coefficients into data vector
        m_data.resize(m_nr_elements * get_nr_coeffs());
        assert(dims[0]*dims[1]*dims[2]*dims[3]*dims[4]==m_data.size());
        H5::DataType data_type = dataset.getDataType();
        assert(data_type.getSize() == sizeof(std::complex<double>));
        dataset.read(m_data.data(), data_type, dataspace);
    } catch (H5::FileIException& e) {
        std::stringstream message;
        message << "Could not load dataset for frequency " << dataset_name << " Mhz";
        throw std::runtime_error(message.str());
    }
}

size_t Dataset::get_index(
    const unsigned int element,
    const unsigned int l,
    const unsigned int m) const
{
    return element * m_nr_pols * m_nr_tetm * m_l_max * m_l_max +
                                                   l * m_l_max +
                                                       m;
}

size_t Dataset::get_nr_coeffs() const
{
    return m_nr_pols * m_nr_tetm * m_l_max * m_l_max;
}

void Dataset::print_alpha(
    const unsigned int element)
{
    const int l_max = get_l_max();

    for (int l = 1; l <= l_max; ++l) {
        for (int abs_m = l; abs_m >=0; --abs_m) {
            auto alpha_ptr = get_alpha_ptr(element, l, abs_m);
            std::cout << *alpha_ptr;
            if (abs_m > 0) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
std::complex<double>* Dataset::get_alpha_ptr(
    const unsigned int element,
    const unsigned int l,
    const unsigned int m)
{
    assert(element < get_nr_elements());
    assert(l <= get_l_max());
    size_t index = get_index(element, l, m);
    return m_data.data() + index;
}

DataFile::DataFile(
    std::string& filename)
{
    // Open file
    std::cout << "read oskar datafile: " << filename << std::endl;
    m_h5_file.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

    // Disable HDF5 error prints
    H5::Exception::dontPrint();
};

std::shared_ptr<Dataset> DataFile::get(
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
