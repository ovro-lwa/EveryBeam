#ifndef EVERYBEAM_CIRCULARSYMMETRIC_VLABEAM_H_
#define EVERYBEAM_CIRCULARSYMMETRIC_VLABEAM_H_

#include <array>
#include <map>
#include <string>

namespace everybeam {
namespace circularsymmetric {
class VLABeam {
 public:
  static std::array<double, 5> GetCoefficients(const std::string& band_name,
                                               double freq);

 private:
  static std::map<int, std::array<double, 5>> GetCoefficients();
  static std::map<char, double> GetFeedConf();
  static char DetermineFeed(double freq, double freqCenter = 0.0);
  static void LimitFreqForBand(char band, double& freq);
};
}  // namespace circularsymmetric
}  // namespace everybeam
#endif  // EVERYBEAM_CIRCULARSYMMETRIC_VLABEAM_H_