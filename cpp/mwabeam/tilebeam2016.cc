#include "tilebeam2016.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>

namespace everybeam::mwabeam {

TileBeam2016::TileBeam2016(const double* delays, bool frequency_interpolation,
                           const std::string& coeff_path)
    : Beam2016Implementation(delays, nullptr, coeff_path),
      frequency_interpolation_(frequency_interpolation),
      mwa_lattitude_(-26.703319),
      mwa_longitude_(116.67081),
      mwa_height_(377.0) {}

void TileBeam2016::ArrayResponse(casacore::MEpoch& time,
                                 casacore::MPosition& array_position, double ra,
                                 double dec, double frequency,
                                 std::complex<double>* gain) {
  casacore::MeasFrame frame(array_position, time);
  const casacore::MDirection::Ref hadec_ref(casacore::MDirection::HADEC, frame);
  const casacore::MDirection::Ref azelgeo_ref(casacore::MDirection::AZELGEO,
                                              frame);
  const casacore::MDirection::Ref j2000_ref(casacore::MDirection::J2000, frame);
  casacore::MPosition wgs = casacore::MPosition::Convert(
      array_position, casacore::MPosition::WGS84)();
  double arr_lattitude =
      wgs.getValue().getLat();  // ant1Pos.getValue().getLat();

  casacore::MDirection::Convert j2000_to_hadec(j2000_ref, hadec_ref),
      j2000_to_azelgeo(j2000_ref, azelgeo_ref);

  ArrayResponse(ra, dec, j2000_ref, j2000_to_hadec, j2000_to_azelgeo,
                arr_lattitude, frequency, gain);
}

void TileBeam2016::ArrayResponse(
    double ra, double dec, const casacore::MDirection::Ref& j2000_ref,
    casacore::MDirection::Convert& j2000_to_hadec,
    casacore::MDirection::Convert& j2000_to_azelgeo, double arr_lattitude,
    double frequency, std::complex<double>* gain) {
  static const casacore::Unit rad_unit("rad");
  casacore::MDirection image_dir(
      casacore::MVDirection(casacore::Quantity(ra, rad_unit),    // RA
                            casacore::Quantity(dec, rad_unit)),  // DEC
      j2000_ref);

  // convert ra, dec to ha
  casacore::MDirection hadec = j2000_to_hadec(image_dir);
  double ha = hadec.getValue().get()[0];
  double sin_lat = std::sin(arr_lattitude), cos_lat = std::cos(arr_lattitude);
  double sin_dec = std::sin(dec), cos_dec = std::cos(dec);
  double cos_ha = std::cos(ha);
  double zenith_distance =
      std::acos(sin_lat * sin_dec + cos_lat * cos_dec * cos_ha);
  casacore::MDirection azel = j2000_to_azelgeo(image_dir);
  double azimuth = azel.getValue().get()[0];
  ArrayResponse(zenith_distance, azimuth, frequency, gain);
}

/**
 * Get the full Jones matrix response of the tile including the dipole
 * reponse and array factor incorporating any mutual coupling effects
 * from the impedance matrix. freq in Hz.
 */
void TileBeam2016::GetTabulatedResponse(double az, double za, double freq,
                                        std::complex<double>* result) {
  // input are radians -> convert to degrees as implementation class expects :
  double az_deg = az * (180.00 / M_PI);
  double za_deg = za * (180.00 / M_PI);
  JonesMatrix jones = CalcJones(az_deg, za_deg, freq, 1);
  result[0] = jones.j00;
  result[1] = jones.j01;
  result[2] = jones.j10;
  result[3] = jones.j11;
}
}  // namespace everybeam::mwabeam
