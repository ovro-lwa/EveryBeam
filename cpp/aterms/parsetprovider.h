// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERMS_PARSETPROVIDER_H_
#define EVERYBEAM_ATERMS_PARSETPROVIDER_H_

#include <string>

namespace everybeam {
namespace aterms {

/**
 * @brief Pure abstract class that can be used to interface the parset (DP3)
 * or configuration file (wsclean) settings for the aterms.
 *
 */
class ParsetProvider {
 public:
  virtual ~ParsetProvider(){};

  virtual std::string GetString(const std::string& key) const = 0;
  virtual std::string GetStringOr(const std::string& key,
                                  const std::string& or_value) const = 0;

  virtual std::vector<std::string> GetStringList(
      const std::string& key) const = 0;

  virtual double GetDoubleOr(const std::string& key, double or_value) const = 0;
  virtual bool GetBool(const std::string& key) const = 0;
  virtual bool GetBoolOr(const std::string& key, bool or_value) const = 0;
};

}  // namespace aterms
}  // namespace everybeam
#endif