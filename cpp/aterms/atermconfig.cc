// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "atermconfig.h"

#include "atermbeam.h"
#include "dldmaterm.h"
#include "everybeamaterm.h"
#include "fitsaterm.h"
#include "pafbeamterm.h"
#include "h5parmaterm.h"

#include "../load.h"
#include "../options.h"
#include "parsetprovider.h"

#include <aocommon/matrix2x2.h>
#include <aocommon/radeccoord.h>

#include <casacore/ms/MeasurementSets/MSAntennaColumns.h>
#include <casacore/tables/Tables/TableKeyword.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <algorithm>

using everybeam::ATermSettings;
using everybeam::coords::CoordinateSystem;

namespace everybeam {
namespace aterms {
void ATermConfig::Read(const casacore::MeasurementSet& ms,
                       const ParsetProvider& reader,
                       const std::string& ms_filename) {
  std::vector<std::string> aterms = reader.GetStringList("aterms");

  if (aterms.empty())
    throw std::runtime_error(
        "No a-term correction given in parset (aterms key is an empty list)");

  for (const std::string& aterm_name : aterms) {
    // Allows to "alias" the aterm type.
    std::string aterm_type =
        reader.GetStringOr(aterm_name + ".type", aterm_name);

    if (aterm_type == "tec") {
      std::vector<std::string> tec_files =
          reader.GetStringList(aterm_name + ".images");
      std::unique_ptr<FitsATerm> f(new FitsATerm(
          n_antennas_, coordinate_system_, settings_.max_support));
      f->OpenTECFiles(tec_files);
      std::string window_str =
          reader.GetStringOr(aterm_name + ".window", "raised-hann");
      aocommon::WindowFunction::Type window =
          aocommon::WindowFunction::GetType(window_str);
      if (window == aocommon::WindowFunction::Tukey) {
        f->SetTukeyWindow(double(settings_.padded_image_width) /
                          settings_.trimmed_image_width);
      } else {
        f->SetWindow(window);
      }
      f->SetDownSample(reader.GetBoolOr(aterm_name + ".downsample", true));
      aterms_.emplace_back(std::move(f));
    } else if (aterm_type == "diagonal") {
      std::vector<std::string> diag_files =
          reader.GetStringList(aterm_name + ".images");
      std::unique_ptr<FitsATerm> f(new FitsATerm(
          n_antennas_, coordinate_system_, settings_.max_support));
      f->OpenDiagGainFiles(diag_files);
      std::string window_str =
          reader.GetStringOr(aterm_name + ".window", "raised-hann");
      aocommon::WindowFunction::Type window =
          aocommon::WindowFunction::GetType(window_str);
      if (window == aocommon::WindowFunction::Tukey) {
        f->SetTukeyWindow(double(settings_.padded_image_width) /
                          settings_.trimmed_image_width);
      } else {
        f->SetWindow(window);
      }
      f->SetDownSample(reader.GetBoolOr(aterm_name + ".downsample", true));
      aterms_.emplace_back(std::move(f));
    } else if (aterm_type == "dldm") {
      std::vector<std::string> dldm_files =
          reader.GetStringList(aterm_name + ".images");
      std::unique_ptr<DLDMATerm> f(new DLDMATerm(
          n_antennas_, coordinate_system_, settings_.max_support));
      f->Open(dldm_files);
      f->SetUpdateInterval(
          reader.GetDoubleOr(aterm_name + ".update_interval", 5.0 * 60.0));
      std::string window_str =
          reader.GetStringOr(aterm_name + ".window", "raised-hann");
      aocommon::WindowFunction::Type window =
          aocommon::WindowFunction::GetType(window_str);
      if (window == aocommon::WindowFunction::Tukey) {
        f->SetTukeyWindow(double(settings_.padded_image_width) /
                          settings_.trimmed_image_width);
      } else {
        f->SetWindow(window);
      }
      f->SetDownSample(reader.GetBoolOr(aterm_name + ".downsample", true));
      aterms_.emplace_back(std::move(f));
    } else if (aterm_type == "beam") {
      bool frequency_interpolation =
          reader.GetBoolOr(aterm_name + ".frequency_interpolation", true);
      bool differential = reader.GetBoolOr(aterm_name + ".differential", false);
      bool use_channel_frequency =
          reader.GetBoolOr(aterm_name + ".usechannelfreq", true);
      std::string element_response_model =
          reader.GetStringOr(aterm_name + ".element_response_model", "default");

      std::unique_ptr<ATermBeam> beam = GetATermBeam(
          ms, coordinate_system_, settings_, frequency_interpolation,
          differential, use_channel_frequency, element_response_model);
      double update_interval = reader.GetDoubleOr(
          aterm_name + ".update_interval", settings_.aterm_update_interval);
      beam->SetUpdateInterval(update_interval);
      aterms_.emplace_back(std::move(beam));
    } else if (aterm_type == "paf") {
      // filenames and ms_filename only needed here.
      if (settings_.filenames.empty() || ms_filename.empty()) {
        throw std::runtime_error(
            "Filenames and ms_filename should be specified for paf aterm "
            "type.");
      }
      auto iter = std::find(settings_.filenames.begin(),
                            settings_.filenames.end(), ms_filename);
      assert(iter != settings_.filenames.end());
      size_t filename_index = iter - settings_.filenames.begin();
      std::vector<std::string> antenna_map =
          reader.GetStringList(aterm_name + ".antenna_map");
      if (antenna_map.size() != n_antennas_) {
        throw std::runtime_error(
            "Antenna map in paf term of aterm config contains " +
            std::to_string(antenna_map.size()) +
            " antennas, whereas the measurement set consists of " +
            std::to_string(n_antennas_));
      }
      std::vector<std::string> beam_map =
          reader.GetStringList(aterm_name + ".beam_map");
      if (beam_map.size() != settings_.filenames.size()) {
        throw std::runtime_error(
            "Number of beams specified in aterm config (" +
            std::to_string(beam_map.size()) +
            ") should match the number of measurement sets specified on the "
            "command line (" +
            std::to_string(settings_.filenames.size()) + ")");
      }
      std::vector<std::string> beam_pointings =
          reader.GetStringList(aterm_name + ".beam_pointings");
      if (beam_pointings.size() != settings_.filenames.size() * 2) {
        throw std::runtime_error("Size of beam pointings is invalid");
      }
      double beam_ra =
          aocommon::RaDecCoord::ParseRA(beam_pointings[filename_index * 2]);
      double beam_dec = aocommon::RaDecCoord::ParseDec(
          beam_pointings[filename_index * 2 + 1]);
      std::string fileTemplate =
          reader.GetString(aterm_name + ".file_template");
      std::unique_ptr<PAFBeamTerm> f(
          new PAFBeamTerm(coordinate_system_, settings_.max_support));
      f->Open(fileTemplate, antenna_map, beam_map[filename_index], beam_ra,
              beam_dec);
      std::string window_str =
          reader.GetStringOr(aterm_name + ".window", "raised-hann");
      aocommon::WindowFunction::Type window =
          aocommon::WindowFunction::GetType(window_str);
      if (window == aocommon::WindowFunction::Tukey) {
        f->SetTukeyWindow(double(settings_.padded_image_width) /
                          settings_.trimmed_image_width);
      } else {
        f->SetWindow(window);
      }
      f->SetDownSample(reader.GetBoolOr(aterm_name + ".downsample", true));
      f->SetReferenceFrequency(
          reader.GetDoubleOr(aterm_name + ".reference_frequency", 0.0));
      aterms_.emplace_back(std::move(f));
    } else if (aterm_type == "h5parm") {
      std::vector<std::string> h5parm_files =
          reader.GetStringList(aterm_name + ".files");

      // Extract antenna names from MS
      std::vector<std::string> station_names(n_antennas_);
      casacore::ScalarColumn<casacore::String> stations(
          ms.antenna(),
          ms.antenna().columnName(casacore::MSAntennaEnums::NAME));

      if (stations.nrow() != n_antennas_) {
        throw std::runtime_error(
            "Number of stations read from measurement set (" +
            std::to_string(stations.nrow()) +
            ") should match the number of stations set in the constructor "
            "command line (" +
            std::to_string(n_antennas_) + ")");
      }

      for (size_t i = 0; i < n_antennas_; ++i) {
        station_names[i] = stations(i);
      }

      std::unique_ptr<H5ParmATerm> f(
          new H5ParmATerm(station_names, coordinate_system_));
      f->Open(h5parm_files);
      f->SetUpdateInterval(reader.GetDoubleOr(aterm_name + ".update_interval",
                                              settings_.aterm_update_interval));
      aterms_.emplace_back(std::move(f));
    } else {
      throw std::runtime_error("The specified aterm type " + aterm_type +
                               " is not recognized");
    }
    aterms_.back()->SetSaveATerms(
        false, settings_.save_aterms_prefix);  // done by config after combining
  }
  if (aterms_.empty()) {
    throw std::runtime_error(
        "The specified a-term configuration does not define any terms to "
        "apply");
  }
  if (aterms_.size() > 1) {
    previous_aterm_values_.resize(aterms_.size());
    for (aocommon::UVector<std::complex<float>>& buf : previous_aterm_values_)
      buf.resize(coordinate_system_.width * coordinate_system_.height *
                 n_antennas_ * 4);
  }
}

bool ATermConfig::Calculate(std::complex<float>* buffer, double time,
                            double frequency, size_t field_id,
                            const double* uvw_in_m) {
  bool is_updated;
  if (aterms_.size() == 1) {
    is_updated =
        aterms_.front()->Calculate(buffer, time, frequency, field_id, uvw_in_m);
  } else {
    is_updated = false;
    for (size_t i = 0; i != aterms_.size(); ++i) {
      is_updated |= aterms_[i]->Calculate(previous_aterm_values_[i].data(),
                                          time, frequency, field_id, uvw_in_m);
    }

    if (is_updated) {
      std::copy(previous_aterm_values_[0].begin(),
                previous_aterm_values_[0].end(), buffer);
      for (size_t i = 1; i != aterms_.size(); ++i) {
        for (size_t j = 0; j != coordinate_system_.width *
                                    coordinate_system_.height * n_antennas_ * 4;
             j += 4) {
          std::complex<float> scratch[4];
          aocommon::Matrix2x2::ATimesB(scratch, &previous_aterm_values_[i][j],
                                       &buffer[j]);
          aocommon::Matrix2x2::Assign(&buffer[j], scratch);
        }
      }
    }
  }

  if (is_updated) {
    SaveATermsIfNecessary(buffer, n_antennas_, coordinate_system_.width,
                          coordinate_system_.height);
  }
  return is_updated;
}

std::unique_ptr<ATermBeam> ATermConfig::GetATermBeam(
    const casacore::MeasurementSet& ms,
    const CoordinateSystem& coordinate_system, const ATermSettings& settings,
    bool frequency_interpolation, bool use_differential_beam,
    bool use_channel_frequency, const std::string& element_response_model) {
  everybeam::Options options = ConvertToEBOptions(
      ms, settings, frequency_interpolation, use_differential_beam,
      use_channel_frequency, element_response_model);
  return std::unique_ptr<ATermBeam>(
      new EveryBeamATerm(ms, coordinate_system, options));
}

everybeam::Options ATermConfig::ConvertToEBOptions(
    const casacore::MeasurementSet& ms, const ATermSettings& settings,
    bool frequency_interpolation, bool use_differential_beam,
    bool use_channel_frequency, const std::string& element_response_model) {
  everybeam::Options options;
  // MWA related
  if (everybeam::GetTelescopeType(ms) ==
      everybeam::TelescopeType::kMWATelescope) {
    // wsclean or dp3 should provide full path to (mwa) coefficients file
    options.coeff_path = settings.coeff_path;
    options.frequency_interpolation = frequency_interpolation;
  }

  everybeam::ElementResponseModel element_response_enum =
      GetElementResponseEnum(element_response_model);

  options.data_column_name = settings.data_column_name;
  options.use_differential_beam = use_differential_beam;
  options.use_channel_frequency = use_channel_frequency;
  options.element_response_model = element_response_enum;
  return options;
}
}  // namespace aterms
}  // namespace everybeam
