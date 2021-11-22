// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MWABEAM_FACTORIALTABLE_H_
#define EVERYBEAM_MWABEAM_FACTORIALTABLE_H_

#include <boost/math/special_functions/factorials.hpp>

#include <aocommon/uvector.h>

namespace everybeam {
namespace mwabeam {
class FactorialTable {
 public:
  FactorialTable(size_t nprecalculated) : table_(nprecalculated) {
    for (size_t i = 0; i != nprecalculated; i++)
      table_[i] = boost::math::factorial<double>(i);
  }

  double operator()(size_t n) const {
    if (n >= table_.size()) {
      return boost::math::factorial<double>(n);
    } else {
      return table_[n];
    }
  }

 private:
  aocommon::UVector<double> table_;
};
}  // namespace mwabeam
}  // namespace everybeam
#endif  // EVERYBEAM_MWABEAM_FACTORIALTABLE_H_
