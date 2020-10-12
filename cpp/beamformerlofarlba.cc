#include "beamformerlofarlba.h"

namespace everybeam {

Antenna::Ptr BeamFormerLofarLBA::Clone() const {
  auto beamformer_clone = std::make_shared<BeamFormerLofarLBA>(
      BeamFormerLofarLBA(coordinate_system_, phase_reference_position_));

  // NOTE: this is an incomplete clone, only creating a deep-copy of the
  // element. In fact, it also hides an upcast from an ElementHamaker into
  // an Element object.
  // The sole and single purpose of Clone() is to be used in
  // Station::SetAntenna!
  Element element_copy = *element_;
  beamformer_clone->SetElement(std::make_shared<Element>(element_copy));
  return beamformer_clone;
}

diag22c_t BeamFormerLofarLBA::LocalArrayFactor(real_t time, real_t freq,
                                               const vector3r_t &direction,
                                               const Options &options) const {
  // Compute the array factor of the field
  diag22c_t array_factor_field = FieldArrayFactor(
      time, freq, direction, options, element_positions_, element_enabled_);

  return array_factor_field;
}
}  // namespace everybeam