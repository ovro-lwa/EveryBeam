// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dldmaterm.h"

#include <aocommon/banddata.h>
#include <aocommon/imagecoordinates.h>

namespace everybeam {
namespace aterms {

DLDMATerm::DLDMATerm(size_t n_antennas,
                     const aocommon::CoordinateSystem& coordinate_system,
                     size_t max_support)
    : FitsATermBase(n_antennas, coordinate_system, max_support),
      update_interval_(60),
      previous_time_(0) {}

void DLDMATerm::Open(const std::vector<std::string>& filenames) {
  readers_.reserve(filenames.size());
  for (const std::string& filename : filenames) {
    readers_.emplace_back(filename, true, true);
    if (readers_.back().NMatrixElements() != 2) {
      throw std::runtime_error(
          "FITS file for dl,dm offsets did not have 2 matrix elements in it");
    }
  }
  InitializeFromFiles(readers_);
}

bool DLDMATerm::Calculate(std::complex<float>* buffer, double time,
                          double frequency, size_t, const double* uvw_in_m) {
  size_t time_index;
  bool requires_recalculation;
  bool position_changed = FindFilePosition(buffer, time, frequency, time_index,
                                           requires_recalculation);
  bool outdated = std::fabs(time - previous_time_) > update_interval_;
  if (!position_changed && !outdated) {
    return false;
  } else {
    if (requires_recalculation || outdated) {
      previous_time_ = time;
      ReadImages(buffer, time_index, frequency, uvw_in_m);
      StoreInCache(frequency, buffer);
    }
    return true;
  }
}

void DLDMATerm::ReadImages(std::complex<float>* buffer, size_t time_index,
                           double frequency, const double* uvw_in_m) {
  const size_t freq_index =
      round((frequency - readers_.front().FrequencyDimensionStart()) /
            readers_.front().FrequencyDimensionIncr());
  const size_t img_index =
      GetTimestep(time_index).img_index * NFrequencies() + freq_index;
  aocommon::FitsReader& reader = readers_[GetTimestep(time_index).reader_index];
  scratch_.resize(GetResampler().ScratchASize());
  dl_image_.resize(GetResampler().ScratchBSize(reader));
  dm_image_.resize(GetResampler().ScratchBSize(reader));
  ReadAndResample(reader, img_index * 2, scratch_, dl_image_, 1.0);
  ReadAndResample(reader, img_index * 2 + 1, scratch_, dm_image_, 1.0);

  double wavel = aocommon::BandData::FrequencyToLambda(frequency);
  for (size_t antenna_idx = 0; antenna_idx != NAntennas(); ++antenna_idx) {
    double uvw[3] = {uvw_in_m[antenna_idx * 3] / wavel,
                     uvw_in_m[antenna_idx * 3 + 1] / wavel,
                     uvw_in_m[antenna_idx * 3 + 2] / wavel};

    std::complex<float>* aterm_buffer =
        buffer + antenna_idx * Width() * Height() * 4;
    EvaluateDLDM(aterm_buffer, dl_image_.data(), dm_image_.data(), uvw);
  }
}

void DLDMATerm::EvaluateDLDM(std::complex<float>* dest, const float* dl,
                             const float* dm, const double* uvw_in_l) {
  using aocommon::ImageCoordinates;

  // For a single source at l,m, we have:
  //   dl = oldl - newl
  //      oldl: measured l
  //      newl: "real" (model) l
  //   V(u,v,w) = I(l, m) exp 2pi i ( lu + mv + nw ), with dn = sqrt(1-l^2-m^2)
  //   - sqrt(1-(l+dl)^2-(m+dm)^2) ~ 0;
  // dl,dm are the shifts in l,m. Given a0 as reference antenna, with u0, v0, w0
  // as coords,
  //   dphase = phase[ I(l, m) exp -2pi i ( dl(u-u0) + dm(v-v0) + dn(w-w0) ) ]
  //          = -2pi i ( l(u-u0) + m(v-v0) + dn(w-w0) )
  // The baselines uvw are already referenced to the first antenna (i.e.
  // uvw=0,0,0 for antenna 0), so uvw_in_l[0] is (u-u0).
  const double u = uvw_in_l[0];
  const double v = uvw_in_l[1];
  const double w = uvw_in_l[2];
  const aocommon::CoordinateSystem& cs = GetCoordinateSystem();
  for (size_t y = 0; y != Height(); ++y) {
    for (size_t x = 0; x != Width(); ++x) {
      double l, m;
      ImageCoordinates::XYToLM(x, y, cs.dl, cs.dm, cs.width, cs.height, l, m);
      l += cs.l_shift;
      m += cs.m_shift;
      const double lproj = l + (*dl);
      const double mproj = m + (*dm);
      const double lm_sq = l * l + m * m;
      const double lmproj_sq = lproj * lproj + mproj * mproj;
      const double dn =
          (lm_sq >= 1.0 || lmproj_sq >= 1.0)
              ? 0.0
              : std::sqrt(1.0 - lmproj_sq) - std::sqrt(1.0 - lm_sq);
      dest[0] = std::polar(1.0, 2.0 * M_PI * (u * (*dl) + v * (*dm) + w * dn));
      dest[1] = 0.0;
      dest[2] = 0.0;
      dest[3] = dest[0];

      ++dl;
      ++dm;
      dest += 4;
    }
  }
}

}  // namespace aterms
}  // namespace everybeam
