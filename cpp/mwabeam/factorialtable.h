#ifndef EVERYBEAM_MWABEAM_FACTORIALTABLE_H_
#define EVERYBEAM_MWABEAM_FACTORIALTABLE_H_

#include <boost/math/special_functions/factorials.hpp>

#include <aocommon/uvector.h>

namespace everybeam::mwabeam {
class FactorialTable {
 public:
  FactorialTable(size_t nprecalculated) : table_(nprecalculated) {
    for (unsigned i = 0; i <= nprecalculated; i++)
      table_[i] = boost::math::factorial<double>(i);
  }

  double operator()(unsigned n) const {
    if (n >= table_.size()) {
      return boost::math::factorial<double>(n);
    }

    if (n < table_.size())
      return table_[n];
    else
      return boost::math::factorial<double>(n);
  }

 private:
  aocommon::UVector<double> table_;
};
}  // namespace everybeam::mwabeam
#endif  // EVERYBEAM_MWABEAM_FACTORIALTABLE_H_
