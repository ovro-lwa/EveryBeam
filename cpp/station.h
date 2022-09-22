// station.h: Representation of the station beam former.
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_STATION_H
#define EVERYBEAM_STATION_H

// \file
// Representation of the station beam former.

#include "elementresponse.h"
#include "antenna.h"
#include "beamformer.h"
#include "coords/itrfdirection.h"
#include "common/types.h"
#include "correctionmode.h"
#include "options.h"
#include <memory>
#include <vector>

#include <aocommon/matrix2x2diag.h>
#include <aocommon/matrix2x2.h>

namespace everybeam {

class Station {
 public:
  /*!
   *  \brief Construct a new Station instance.
   *
   *  \param name Name of the station.
   *  \param position Position of the station (ITRF, m).
   */
  Station(const std::string& name, const vector3r_t& position,
          const Options& options = Options());

  //! Return the name of the station.
  const std::string& GetName() const;

  //! Return the position of the station (ITRF, m).
  const vector3r_t& GetPosition() const;

  /*!
   *  \brief Set the phase reference position. This is the position where the
   *  delay of the incoming plane wave is assumed to be zero.
   *
   *  \param reference Phase reference position (ITRF, m).
   *
   *  By default, it is assumed the position of the station is also the phase
   *  reference position. Use this method to set the phase reference position
   *  explicitly when this assumption is false.
   */
  void SetPhaseReference(const vector3r_t& reference);

  //! Return the phase reference position (ITRF, m). \see
  //! Station::setPhaseReference()
  const vector3r_t& GetPhaseReference() const;

  /**
   *  \brief Compute the full response of the station for a plane wave of
   * frequency \p freq, arriving from direction \p direction, with the %station
   * beam former steered towards \p station0, and, for HBA stations, the analog
   *  %tile beam former steered towards \p tile0. For LBA stations, \p tile0
   *  has no effect.
   *
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Frequency of the plane wave (Hz).
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 %Station beam former reference frequency (Hz).
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \return Jones matrix that represents the %station response.
   *
   *  For any given sub-band, the (%LOFAR) station beam former computes weights
   *  for a single reference frequency. Usually, this reference frequency is
   *  the center frequency of the sub-band. For any frequency except the
   *  reference frequency, these weights are an approximation. This aspect of
   *  the system is taken into account in the computation of the response.
   *  Therefore, both the frequency of interest \p freq and the reference
   *  frequency \p freq0 need to be specified.
   *
   *  The directions \p direction, \p station0, and \p tile0 are vectors that
   *  represent a direction of \e arrival. These vectors have unit length and
   *  point \e from the ground \e towards the direction from which the plane
   *  wave arrives.
   */
  aocommon::MC2x2 Response(real_t time, real_t freq,
                           const vector3r_t& direction, real_t freq0,
                           const vector3r_t& station0, const vector3r_t& tile0,
                           const bool rotate = true) const;

  /**
   * \brief This method is similar to the above Response() function, but adds
   * parameter \p mode to request the response for a specific mode.
   * \see CorrectionMode.
   */
  aocommon::MC2x2 Response(CorrectionMode mode, real_t time, real_t freq,
                           const vector3r_t& direction, real_t freq0,
                           const vector3r_t& station0, const vector3r_t& tile0,
                           const bool rotate = true) const {
    switch (mode) {
      case CorrectionMode::kNone:
        return aocommon::MC2x2::Unity();
      case CorrectionMode::kFull:
        return Response(time, freq, direction, freq0, station0, tile0, rotate);
      case CorrectionMode::kArrayFactor:
        return ArrayFactor(time, freq, direction, freq0, station0, tile0);
      case CorrectionMode::kElement:
        return ComputeElementResponse(time, freq, direction, false, rotate);
    }
    throw std::runtime_error("Invalid mode");
  }

  /*!
   *  \brief Compute the array factor of the station for a plane wave of
   *  frequency \p freq, arriving from direction \p direction, with the
   *  %station beam former steered towards \p station0, and, for HBA stations
   *  the analog %tile beam former steered towards \p tile0. For LBA stations,
   *  \p tile0 has no effect.
   *
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Frequency of the plane wave (Hz).
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 %Station beam former reference frequency (Hz).
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \return A diagonal matrix with the array factor of the X and Y antennae.
   *
   *  For any given sub-band, the (%LOFAR) station beam former computes weights
   *  for a single reference frequency. Usually, this reference frequency is
   *  the center frequency of the sub-band. For any frequency except the
   *  reference frequency, these weights are an approximation. This aspect of
   *  the system is taken into account in the computation of the response.
   *  Therefore, both the frequency of interest \p freq and the reference
   *  frequency \p freq0 need to be specified.
   *
   *  The directions \p direction, \p station0, and \p tile0 are vectors that
   *  represent a direction of \e arrival. These vectors have unit length and
   *  point \e from the ground \e towards the direction from which the plane
   *  wave arrives.
   */
  aocommon::MC2x2Diag ArrayFactor(real_t time, real_t freq,
                                  const vector3r_t& direction, real_t freq0,
                                  const vector3r_t& station0,
                                  const vector3r_t& tile0) const;

  /*!
   *  \name Convenience member functions
   *  These member functions perform the same function as the corresponding
   *  non-template member functions, for a list of frequencies or (frequency,
   *  reference frequency) pairs.
   */
  // @{

  /*!
   *  \brief Convenience method to compute the response of the station for a
   *  list of frequencies, and a fixed reference frequency.
   *
   *  \param count Number of frequencies.
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Input iterator for a list of frequencies (Hz) of length
   *  \p count.
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 %Station beam former reference frequency (Hz).
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \param buffer Output iterator with room for \p count instances of type
   *  ::aocommon::MC2x2.
   *
   *  \see response(real_t time, real_t freq, const vector3r_t &direction,
   *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
   */
  template <typename T, typename U>
  void Response(unsigned int count, real_t time, T freq,
                const vector3r_t& direction, real_t freq0,
                const vector3r_t& station0, const vector3r_t& tile0, U buffer,
                const bool rotate = true) const;

  /*!
   *  \brief Convenience method to compute the array factor of the station for
   *  a list of frequencies, and a fixed reference frequency.
   *
   *  \param count Number of frequencies.
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Input iterator for a list of frequencies (Hz) of length
   *  \p count.
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 %Station beam former reference frequency (Hz).
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \param buffer Output iterator with room for \p count instances of type
   *  ::aocommon::MC2x2.
   *
   *  \see ArrayFactor(real_t time, real_t freq, const vector3r_t &direction,
   *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
   */
  template <typename T, typename U>
  void ArrayFactor(unsigned int count, real_t time, T freq,
                   const vector3r_t& direction, real_t freq0,
                   const vector3r_t& station0, const vector3r_t& tile0,
                   U buffer) const;

  /*!
   *  \brief Convenience method to compute the response of the station for a
   *  list of (frequency, reference frequency) pairs.
   *
   *  \param count Number of frequencies.
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Input iterator for a list of frequencies (Hz) of length
   *  \p count.
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 Input iterator for a list of %Station beam former reference
   *  frequencies (Hz) of length \p count.
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \param buffer Output iterator with room for \p count instances of type
   *  ::aocommon::MC2x2.
   *
   *  \see response(real_t time, real_t freq, const vector3r_t &direction,
   *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
   */
  template <typename T, typename U>
  void Response(unsigned int count, real_t time, T freq,
                const vector3r_t& direction, T freq0,
                const vector3r_t& station0, const vector3r_t& tile0, U buffer,
                const bool rotate = true) const;

  /*!
   *  \brief Convenience method to compute the array factor of the station for
   *  list of (frequency, reference frequency) pairs.
   *
   *  \param count Number of frequencies.
   *  \param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   *  \param freq Input iterator for a list of frequencies (Hz) of length
   *  \p count.
   *  \param direction Direction of arrival (ITRF, m).
   *  \param freq0 %Station beam former reference frequency (Hz).
   *  \param station0 %Station beam former reference direction (ITRF, m).
   *  \param tile0 Tile beam former reference direction (ITRF, m).
   *  \param rotate Boolean deciding if paralactic rotation should be applied.
   *  \param buffer Output iterator with room for \p count instances of type
   *  ::aocommon::MC2x2.
   *
   *  \see ArrayFactor(real_t time, real_t freq, const vector3r_t &direction,
   *  real_t freq0, const vector3r_t &station0, const vector3r_t &tile0) const
   */
  template <typename T, typename U>
  void ArrayFactor(unsigned int count, real_t time, T freq,
                   const vector3r_t& direction, T freq0,
                   const vector3r_t& station0, const vector3r_t& tile0,
                   U buffer) const;

  // @}

  //! Returns a pointer to the ElementResponse class
  std::shared_ptr<const ElementResponse> GetElementResponse() const {
    return element_response_;
  }

  /**
   * @brief Compute the Jones matrix for the element response
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival. If is_local is true: (ENU, m) else
   * direction vector in global coord system is assumed.
   * @param is_local Use local east-north-up system (true) or global coordinate
   * system (false).
   * @param id Element id
   * @param rotate Boolean deciding if paralactic rotation should be applied.
   * @return aocommon::MC2x2 Jones matrix of element response
   */
  aocommon::MC2x2 ComputeElementResponse(real_t time, real_t freq,
                                         const vector3r_t& direction, size_t id,
                                         bool is_local, bool rotate) const;

  /**
   * @brief Compute the Jones matrix for the element response
   *
   * @param time Time, modified Julian date, UTC, in seconds (MJD(UTC), s).
   * @param freq Frequency of the plane wave (Hz).
   * @param direction Direction of arrival. If is_local is true: (ENU, m) else
   * direction vector in global coord system is assumed.
   * @param is_local Use local east-north-up system (true) or global coordinate
   * system (false).
   * @param rotate Boolean deciding if paralactic rotation should be applied.
   * @return aocommon::MC2x2 Jones matrix of element response
   */
  aocommon::MC2x2 ComputeElementResponse(real_t time, real_t freq,
                                         const vector3r_t& direction,
                                         bool is_local, bool rotate) const;

  //! Specialized implementation of response function.
  aocommon::MC2x2 Response(real_t time, real_t freq,
                           const vector3r_t& direction) const {
    return antenna_->Response(*element_response_, time, freq, direction);
  }

  //! Set antenna attribute, usually a BeamFormer, but can also be an Element
  void SetAntenna(std::shared_ptr<Antenna> antenna);

  std::shared_ptr<Antenna> GetAntenna() const { return antenna_; }

 private:
  vector3r_t NCP(real_t time) const;
  vector3r_t NCPPol0(real_t time) const;

  std::string name_;
  vector3r_t position_;
  Options options_;
  vector3r_t phase_reference_;
  std::shared_ptr<const ElementResponse> element_response_;
  // element_ either refers to antenna_ or an Element inside antenna_.
  // Besides Station, no one has (shared) ownership of antenna_.
  std::shared_ptr<Element> element_;
  std::shared_ptr<Antenna> antenna_;

  coords::ITRFDirection::Ptr ncp_;
  /** Reference direction for NCP observations.
   *
   * NCP pol0 is the direction used as reference in the coordinate system
   * when the target direction is close to/at the NCP. The regular coordinate
   * system rotates local east to that defined with respect to the NCP,
   * which is undefined at the NCP.
   * It is currently defined as ITRF position (1.0, 0.0, 0.0).
   *
   * Added by Maaijke Mevius, December 2018.
   */
  coords::ITRFDirection::Ptr ncp_pol0_;
};

// ------------------------------------------------------------------------- //
// - Implementation: Station                                               - //
// ------------------------------------------------------------------------- //

template <typename T, typename U>
void Station::Response(unsigned int count, real_t time, T freq,
                       const vector3r_t& direction, real_t freq0,
                       const vector3r_t& station0, const vector3r_t& tile0,
                       U buffer, const bool rotate) const {
  for (unsigned int i = 0; i < count; ++i) {
    *buffer++ =
        Response(time, *freq++, direction, freq0, station0, tile0, rotate);
  }
}

template <typename T, typename U>
void Station::ArrayFactor(unsigned int count, real_t time, T freq,
                          const vector3r_t& direction, real_t freq0,
                          const vector3r_t& station0, const vector3r_t& tile0,
                          U buffer) const {
  for (unsigned int i = 0; i < count; ++i) {
    *buffer++ = ArrayFactor(time, *freq++, direction, freq0, station0, tile0);
  }
}

template <typename T, typename U>
void Station::Response(unsigned int count, real_t time, T freq,
                       const vector3r_t& direction, T freq0,
                       const vector3r_t& station0, const vector3r_t& tile0,
                       U buffer, const bool rotate) const {
  for (unsigned int i = 0; i < count; ++i) {
    *buffer++ =
        Response(time, *freq++, direction, *freq0++, station0, tile0, rotate);
  }
}

template <typename T, typename U>
void Station::ArrayFactor(unsigned int count, real_t time, T freq,
                          const vector3r_t& direction, T freq0,
                          const vector3r_t& station0, const vector3r_t& tile0,
                          U buffer) const {
  for (unsigned int i = 0; i < count; ++i) {
    *buffer++ =
        ArrayFactor(time, *freq++, direction, *freq0++, station0, tile0);
  }
}
}  // namespace everybeam
#endif
