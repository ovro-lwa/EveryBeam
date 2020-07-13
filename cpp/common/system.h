#ifndef EVERYBEAM_COMMON_SYSTEM_H_
#define EVERYBEAM_COMMON_SYSTEM_H_

#include <string>
#include <sched.h>

namespace everybeam {
namespace common {

class System {
 public:
  static std::size_t ProcessorCount() {
#ifdef __APPLE__
    return sysconf(_SC_NPROCESSORS_ONLN);
#else
    cpu_set_t cs;
    CPU_ZERO(&cs);
    sched_getaffinity(0, sizeof cs, &cs);

    int count = 0;
    for (int i = 0; i < CPU_SETSIZE; i++) {
      if (CPU_ISSET(i, &cs)) count++;
    }

    return count;
#endif
  }
};

}  // namespace common
}  // namespace everybeam

#endif  // EVERYBEAM_COMMON_SYSTEM_H_