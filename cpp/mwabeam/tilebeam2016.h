// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MWABEAM_TILEBEAM2016_H_
#define EVERYBEAM_MWABEAM_TILEBEAM2016_H_

#include <complex>
#include <map>
#include <set>

#include "beam2016implementation.h"

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>

namespace everybeam {
namespace mwabeam {

class TileBeam2016 : public Beam2016Implementation {
 public:
  TileBeam2016(const double* delays, bool frequency_interpolation,
               const std::string& coeff_path);

  /**
   * @brief API method for computing MWA array response
   *
   * @param ra right ascension (rad)
   * @param dec declination (rad)
   * @param j2000_ref J2000 ref coordinates
   * @param j2000_to_hadecref HADEC coordinates
   * @param j2000_to_azelgeoref AZELGEO coordinates
   * @param arr_lattitude Lattitude
   * @param frequency Frequency (Hz)
   * @param gain Gain matrix
   */
  void ArrayResponse(double ra, double dec,
                     const casacore::MDirection::Ref& j2000_ref,
                     casacore::MDirection::Convert& j2000_to_hadecref,
                     casacore::MDirection::Convert& j2000_to_azelgeoref,
                     double arr_lattitude, double frequency,
                     std::complex<double>* gain);

  /**
   * @brief Compute MWA array response in given zenith/azimuth direction
   *
   * @param zenith_angle Zenith angle (rad)
   * @param azimuth Azimuthal angle (rad)
   * @param frequency Frequency (Hz)
   * @param gain Gain matrix
   */
  void ArrayResponse(double zenith_angle, double azimuth, double frequency,
                     std::complex<double>* gain) {
    // As yet, this conditional is effectively redundant
    if (frequency_interpolation_)
      GetInterpolatedResponse(azimuth, zenith_angle, frequency, gain);
    else
      GetTabulatedResponse(azimuth, zenith_angle, frequency, gain);
  }

 private:
  bool frequency_interpolation_;

  /**
   * Get the full Jones matrix response of the tile including the dipole
   * reponse and array factor incorporating any mutual coupling effects
   * from the impedance matrix. freq in Hz.
   */
  void GetTabulatedResponse(double az, double za, double freq,
                            std::complex<double>* result);

  /**
   * Create a few tabulated responses and interpolated over these.
   */
  void GetInterpolatedResponse(double az, double za, double freq,
                               std::complex<double>* result) {
    // Not implemented yet: just call normal function
    GetTabulatedResponse(az, za, freq, result);
  }
};
}  // namespace mwabeam
}  // namespace everybeam
#endif  // EVERYBEAM_MWABEAM_TILEBEAM2016_H_
