#include "OSKARElementResponse.h"
#include "oskar.h"
#include "config.h"
#include <iostream>

namespace LOFAR {
namespace StationResponse{

void OSKARElementResponseDipole::response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    double dipole_length_m = 1; // TODO
    std::complex<double>* response_ptr = (std::complex<double> *) response;

    double phi_x = phi;
    double phi_y = phi + M_PI/2;
    oskar_evaluate_dipole_pattern_double(1, &theta, &phi_x, freq, dipole_length_m, response_ptr);
    oskar_evaluate_dipole_pattern_double(1, &theta, &phi_y, freq, dipole_length_m, response_ptr + 2);
}

OSKARElementResponseSphericalWave::OSKARElementResponseSphericalWave()
{
    std::string path = get_path("oskar.h5");
    m_datafile.reset(new DataFile(path));
}

void OSKARElementResponseSphericalWave::response(
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    // This ElementResponse model is element specific, so an element_id is required
    // to know for what element the response needs to be evaluated
    // A std::invalid_argument exception is thrown although strictly speaking
    // it are not the given arguments that are invalid, but the response(...) method with
    // a different signature should have been called.
    throw std::invalid_argument("OSKARElementResponseSphericalWave: missing argument element_id");
}

void OSKARElementResponseSphericalWave::response(
    int element_id,
    double freq,
    double theta,
    double phi,
    std::complex<double> (&response)[2][2]) const
{
    auto dataset = m_datafile->get(freq);
    auto l_max = dataset->get_l_max();

    std::complex<double>* response_ptr = (std::complex<double> *) response;
    std::complex<double>* alpha_ptr = dataset->get_alpha_ptr(element_id);

    double phi_x = phi;
    double phi_y = phi + M_PI/2;
    oskar_evaluate_spherical_wave_sum_double(1, &theta, &phi_x, &phi_y, l_max, alpha_ptr, response_ptr);
}

std::string OSKARElementResponseSphericalWave::get_path(
    const char* filename) const
{
    std::stringstream ss;
    ss << LOFARBEAM_DATA_DIR << "/";
    ss << filename;
    return ss.str();
}

} // namespace StationResponse
} // namespace LOFAR
