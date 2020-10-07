#ifndef EVERYBEAM_BEAMFORMERLOFARHBA_H
#define EVERYBEAM_BEAMFORMERLOFARHBA_H

#include <complex>
#include <vector>
#include <memory>

#include "beamformer.h"
#include "elementhamaker.h"
#include "common/types.h"
#include "antenna.h"

namespace everybeam {
/**
 * @brief Optimized implementation of the BeamFormer class for the LOFAR HBA
 * telescope in combination with Hamaker element response model.
 *
 */
class BeamFormerLofarHBA : public Antenna {
 public:
  /**
   * @brief Construct a new BeamFormer object
   *
   */
  BeamFormerLofarHBA() : Antenna() {}

  /**
   * @brief Construct a new BeamFormerLofarHBA object given a coordinate system.
   *
   * @param coordinate_system
   */
  BeamFormerLofarHBA(const CoordinateSystem &coordinate_system)
      : Antenna(coordinate_system) {}

  /**
   * @brief Construct a new BeamFormerLofarHBA object given a coordinate system
   * and a phase reference position
   *
   * @param coordinate_system
   * @param phase_reference_position
   */
  BeamFormerLofarHBA(CoordinateSystem coordinate_system,
                     vector3r_t phase_reference_position)
      : Antenna(coordinate_system, phase_reference_position) {}

  BeamFormerLofarHBA(vector3r_t phase_reference_position)
      : Antenna(phase_reference_position) {}

  /**
   * @brief Returns an (incomplete!) clone of the BeamFormerLofarHBA class
   * only the element_ is copied. This method is intended to be exclusively
   * used in Station::SetAntenna!
   *
   * @return Antenna::Ptr
   */
  Antenna::Ptr Clone() const override;

  /**
   * @brief Set the (unique) Tile for the BeamFormerLofarHBA object.
   *
   * @param beamformer
   */
  void SetTile(BeamFormer::Ptr beamformer) { tile_ = beamformer; }

  /**
   * @brief Set the (unique) Element object for the BeamFormerLofarHBA object.
   *
   * @param element
   */
  void SetElement(Element::Ptr element) { element_ = element; }
  // void SetElement(Element::Ptr element) { element_ = element; }

  /**
   * @brief Add tile position to the tile_positions_ array
   *
   * @param position
   */
  void AddTilePosition(const vector3r_t &position) {
    tile_positions_.push_back(position);
  }

  /**
   * @brief Mark whether tile is enabled by pushing back boolean array to
   * tile_enabled_ array
   *
   * @param enabled
   */
  void AddTileEnabled(const std::array<bool, 2> enabled) {
    tile_enabled_.push_back(enabled);
  }

  /**
   * @brief Add element position to the element_positions_ array
   *
   * @param position
   */
  void AddElementPosition(const vector3r_t &position) {
    element_positions_.push_back(position);
  }

  /**
   * @brief Returns the element_ object of the BeamFormerLofarHBA
   *
   * @return std::shared_ptr<Element>
   */
  std::shared_ptr<Element> GetElement() const { return element_; };

 private:
  // Compute the geometric response either for the tiles, or for the elements
  // within the BeamFormerLofarHBA. Method assumes that the direction is
  // specified as the (frequency weighted) difference between the pointing_dir
  // and the probing direction
  std::vector<std::complex<double>> ComputeGeometricResponse(
      std::vector<vector3r_t> phase_reference_positions,
      const vector3r_t &direction) const;

  // Array factor
  virtual diag22c_t LocalArrayFactor(real_t time, real_t freq,
                                     const vector3r_t &direction,
                                     const Options &options) const override;

  // Array factor at Field level, weighting the tiles
  virtual diag22c_t FieldArrayFactor(real_t time, real_t freq,
                                     const vector3r_t &direction,
                                     const Options &options) const;

  // Array Factor at tile level, weighting the elements
  virtual std::complex<double> TileArrayFactor(real_t time, real_t freq,
                                               const vector3r_t &direction,
                                               const Options &options) const;

  // Override of the LocalResponse method
  virtual matrix22c_t LocalResponse(real_t time, real_t freq,
                                    const vector3r_t &direction,
                                    const Options &options) const override;

  // BeamFormerLofarHBA has only one unique tile, and ntiles
  // unique tile positions
  std::shared_ptr<BeamFormer> tile_;
  std::vector<vector3r_t> tile_positions_;

  // Is tile enabled?
  std::vector<std::array<bool, 2>> tile_enabled_;

  // Each BeamFormerLofarHBA stores one unique element, and 16 unique positions
  // Usually, the Element will be of type ElementHamaker.
  std::shared_ptr<Element> element_;
  std::vector<vector3r_t> element_positions_;
};
}  // namespace everybeam
#endif
