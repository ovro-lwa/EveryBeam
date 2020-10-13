#ifndef EVERYBEAM_SINGLETON_H
#define EVERYBEAM_SINGLETON_H

namespace everybeam {
namespace common {
template <typename T>
struct Singleton {
  // Factory function to obtain the one-and-only instance of the Singleton
  static std::shared_ptr<T> GetInstance() {
    // Static variable, initialized on first call to GetInstance()
    static std::shared_ptr<T> instance = std::make_shared<T>();
    return instance;
  }

 private:
  // Make the constructor private, to prevent direct instantiation
  // outside of the factory function GetInstance()
  Singleton() {}

  // Forbid to make copies of the Singleton, by deleting the
  // copy and assignment constructors
  Singleton(Singleton const&) = delete;
  void operator=(Singleton const&) = delete;
};
}  // namespace common
}  // namespace everybeam
#endif