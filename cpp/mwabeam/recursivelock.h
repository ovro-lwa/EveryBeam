#ifndef EVERYBEAM_MWABEAM_RECURSIVELOCK_H_
#define EVERYBEAM_MWABEAM_RECURSIVELOCK_H_

#include <cstring>
#include <mutex>
#include <system_error>

namespace everybeam::mwabeam {

template <typename Mutex>
class RecursiveLock {
 public:
  RecursiveLock() : mutex_(nullptr), nlocks_(0) {}

  RecursiveLock(RecursiveLock&& other)
      : mutex_(other.mutex_), nlocks_(other.nlocks_) {
    other.mutex_ = nullptr;
    other.nlocks_ = 0;
  }

  RecursiveLock(Mutex& mutex) : mutex_(&mutex), nlocks_(1) { mutex_->lock(); }

  RecursiveLock(Mutex& mutex, std::defer_lock_t) noexcept
      : mutex_(&mutex), nlocks_(0) {}

  RecursiveLock(Mutex& mutex, std::try_to_lock_t) : mutex_(&mutex), nlocks_(0) {
    try_lock();
  }

  RecursiveLock(Mutex& mutex, std::adopt_lock_t) noexcept
      : mutex_(&mutex), nlocks_(1) {}

  ~RecursiveLock() {
    if (nlocks_ != 0) mutex_->unlock();
  }

  RecursiveLock& operator=(const RecursiveLock& other) {
    if (nlocks_ != 0) mutex_->unlock();
    mutex_ = other.mutex_;
    nlocks_ = other.nLocks;
    other.mutex_ = nullptr;
    other.nlocks_ = 0;
  }

  void lock() {
    if (mutex_ == nullptr)
      throw std::system_error(
          std::make_error_code(std::errc::operation_not_permitted),
          "lock() called on RecursiveLock without mutex");
    if (nlocks_ == 0) mutex_->lock();
    ++nlocks_;
  }

  bool try_lock() {
    if (mutex_ == nullptr)
      throw std::system_error(
          std::make_error_code(std::errc::operation_not_permitted),
          "lock() called on RecursiveLock without mutex");
    if (nlocks_ == 0) {
      if (mutex_->try_lock()) {
        ++nlocks_;
        return true;
      } else {
        return false;
      }
    } else {
      ++nlocks_;
      return true;
    }
  }

  void unlock() {
    nlocks_--;
    if (mutex_ == nullptr)
      throw std::system_error(
          std::make_error_code(std::errc::operation_not_permitted),
          "unlock() called on RecursiveLock without mutex");
    else if (nlocks_ == 0)
      mutex_->unlock();
  }

  void swap(RecursiveLock& other) noexcept {
    std::swap(mutex_, other.mutex_);
    std::swap(nlocks_, other.nlocks_);
  }

  bool owns_lock() const noexcept { return nlocks_ != 0; }

  Mutex* release() noexcept {
    Mutex* m = mutex_;
    mutex_ = nullptr;
    nlocks_ = 0;
    return m;
  }

  Mutex* mutex() const noexcept { return mutex_; }

  operator bool() const noexcept { return owns_lock(); }

 private:
  Mutex* mutex_;
  size_t nlocks_;
};
}  // namespace everybeam::mwabeam
#endif  // EVERYBEAM_MWABEAM_RECURSIVELOCK_H_
