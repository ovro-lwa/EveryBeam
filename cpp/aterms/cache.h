// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_ATERM_CACHE_H
#define EVERYBEAM_ATERM_CACHE_H

#include <algorithm>
#include <complex>
#include <memory>
#include <vector>
#include <limits>

namespace everybeam {
namespace aterms {
/**
 * A simple a-term cache that stores aterms for all frequencies in a single
 * timestep. Since each timestep will have the same frequency values, the cache
 * maintains buffers for each frequency value, and once all frequencies have
 * been stored at least once, no more allocations are performed.
 *
 * Used by the FitsATerm class.
 */
class Cache {
 public:
  /**
   * Construct the cache.
   */
  Cache() : aterm_size_(0) {}
  Cache(size_t aterm_size) : aterm_size_(aterm_size) {}

  Cache(const Cache&) = delete;
  Cache(Cache&&) = default;
  Cache& operator=(const Cache&) = delete;
  Cache& operator=(Cache&&) = default;

  [[gnu::visibility("default")]] static const size_t
      kNotFound;  // = std::numeric_limits<size_t>::max();

  /**
   * Clears all stored values, such that e.g. the cache is ready for a new
   * timestep. After this call, @ref Find() will always return @ref kNotFound .
   */
  void Reset() {
    // The cache is 'emptied', but we don't entirely clear the data vector at
    // the top level, as we can reuse the space from the individual entries
    // and prevent reallocation.
    for (Entry& entry : entries_) {
      entry.isValid = false;
    }
  }

  /**
   * Returns the index of the frequency, or @ref kNotFound if not found.
   */
  size_t Find(double frequency) const {
    auto iter =
        std::lower_bound(frequencies_.begin(), frequencies_.end(), frequency);
    if (iter == frequencies_.end() || *iter != frequency)
      return kNotFound;
    else {
      size_t index = iter - frequencies_.begin();
      if (entries_[index].isValid)
        return index;
      else
        return kNotFound;
    }
  }

  /**
   * Retrieve a buffer given a frequency index (as returned by @ref Find() ).
   * @param destination Array with @ref ATermSize() elements.
   */
  void Get(size_t index, std::complex<float>* destination) const {
    std::copy_n(entries_[index].ptr.get(), aterm_size_, destination);
  }

  /**
   * Store the data for the given frequency in the cache. The data array should
   * be an array with @ref ATermSize() elements.
   */
  void Store(double frequency, const std::complex<float>* data) {
    auto iter =
        std::lower_bound(frequencies_.begin(), frequencies_.end(), frequency);
    if (iter != frequencies_.end() && *iter == frequency) {
      size_t index = iter - frequencies_.begin();
      std::copy_n(data, aterm_size_, entries_[index].ptr.get());
      entries_[index].isValid = true;
    } else {
      size_t index = iter - frequencies_.begin();
      frequencies_.emplace(iter, frequency);
      Entry newEntry;
      newEntry.ptr.reset(new std::complex<float>[aterm_size_]);
      newEntry.isValid = true;
      std::copy_n(data, aterm_size_, newEntry.ptr.get());
      entries_.emplace(entries_.begin() + index, std::move(newEntry));
    }
  }

  /**
   * Size of one aterm buffer (in number of complex float values).
   */
  size_t ATermSize() const { return aterm_size_; }

 private:
  // This array is always kept sorted
  std::vector<double> frequencies_;
  size_t aterm_size_;

  struct Entry {
    std::unique_ptr<std::complex<float>[]> ptr;
    bool isValid;
  };

  // entries_[index] corresponds with frequencies_[index]
  std::vector<Entry> entries_;
};
}  // namespace aterms
}  // namespace everybeam
#endif
