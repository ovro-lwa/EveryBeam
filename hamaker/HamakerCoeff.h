#ifndef HAMAKER_COEFF_H
#define HAMAKER_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

//! Hamaker coefficients
class HamakerCoefficients {
    public:
        //! Default constructor
        HamakerCoefficients();

        //! Constructor for reading coeff from file
        HamakerCoefficients(
            std::string& filename);

        //! Constructor for writing coeff to file
        HamakerCoefficients(
            const double freq_center,
            const double freq_range,
            const unsigned int nHarmonics,
            const unsigned int nPowerTheta,
            const unsigned int nPowerFreq);

        /**
         * @brief Set Hamaker coefficients
         * 
         * @param n 
         * @param t 
         * @param f 
         * @param value 
         */
        void set_coeff(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f,
            std::pair<std::complex<double>, std::complex<double>> value);

        void set_coeffs(
            const std::complex<double>* coeff);

        void set_coeffs(
            const std::vector<std::complex<double>> coeff);

        // Get
        size_t get_nr_coeffs() const;

        double get_freq_center() const { return m_freq_center; }

        double get_freq_range() const { return m_freq_range; }

        unsigned int get_nHarmonics() const { return m_nHarmonics; }

        unsigned int get_nPowerTheta() const { return m_nPowerTheta; }

        unsigned int get_nPowerFreq() const { return m_nPowerFreq; }


        std::pair<std::complex<double>, std::complex<double>> get_coeff(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f);

        // HDF5 I/O
        void read_coeffs(
            std::string& filename);

        void write_coeffs(
            std::string& filename);

        // Debugging
        void print_coeffs();

    private:
        // Methods
        size_t get_index(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f);

        // Parameters
        double m_freq_center;
        double m_freq_range;
        unsigned int m_nHarmonics;
        unsigned int m_nPowerTheta;
        unsigned int m_nPowerFreq;
        const unsigned int m_nInner = 2;

        // Data
        std::vector<std::complex<double>> m_coeff;

        // HDF5
        std::string m_dataset_name = "coeff";
        const unsigned int m_dataset_rank = 4;
};

#endif
