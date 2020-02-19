#ifndef OSKAR_SPHERICAL_WAVE_COEFF_H
#define OSKAR_SPHERICAL_WAVE_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>
#include <memory>
#include <mutex>
#include <map>

#include <H5Cpp.h>

class Dataset {
    public:
        Dataset(
            H5::H5File& h5_file,
            const unsigned int freq);

        // Get
        size_t get_nr_elements() const { return m_nr_elements; };
        size_t get_l_max() const { return m_l_max; };
        size_t get_nr_coeffs() const;

        std::complex<double>* get_alpha_ptr(
            const unsigned int element);

        // Debugging
        void print_alpha(
            const unsigned int element = 0);

    private:
        // Methods
        size_t get_index(
            const unsigned int element,
            const unsigned int l,
            const unsigned int m) const;

        // Constants
        const unsigned int m_dataset_rank = 5;
        const unsigned int m_nr_pols = 2;
        const unsigned int m_nr_tetm = 2;

        // Members
        std::vector<std::complex<double>> m_data;
        unsigned int m_nr_elements;
        unsigned int m_l_max;
};

class DataFile {
    public:
        // Constructor for reading coeff from file
        DataFile(
            std::string& filename);

        std::shared_ptr<Dataset> get(
            const unsigned int freq);

    private:
        // Coeffs;
        std::map<unsigned int, std::shared_ptr<Dataset>> m_map;

        // HDF5
        std::string m_filename;
        std::unique_ptr<H5::H5File> m_h5_file;
        mutable std::mutex m_mutex;
};

#endif
