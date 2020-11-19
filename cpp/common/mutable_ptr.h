// mutableptr.h: Representation of an LBA antenna field.
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVERYBEAM_MUTABLEPTR_H
#define EVERYBEAM_MUTABLEPTR_H

#include <memory>

namespace everybeam {
namespace common {
/*!
 * \brief MutablePtr is a mutable smart pointer derived from std::shared_ptr.
 *
 * Its purpose is to have a pointer that
 *     1. can by copied/passed by value,
 *     2. the original can be modified after it has been copied
 *     3. the copies will dereference to the new value
 *
 * This behaviour is different from a `reset` of a std::shared_ptr
 * There the copies will keep pointing to the original object
 *
 * The inner pointer is a shared_ptr which manages the object it points to
 *
 * Usage:
 *
 * \code
 * #include<iostream>
 * #include "mutable_ptr.h"
 *
 * class Msg {
 * public:
 *     Msg(const std::string &t) : text(t) {}
 *     std::string text;
 * };
 *
 * typedef MutablePtr<Msg> MsgMutPtr;
 *
 * int main()
 * {
 * auto m1 = std::make_shared<Msg>("hoi");
 * auto m2 = std::make_shared<Msg>("hallo");
 *
 * MsgPtr mut_ptr(m1);
 * auto mut_ptr_copy = mut_ptr;
 *
 * mut_ptr.set(m2);
 *
 * std::cout << mut_ptr->text << std::endl;
 * std::cout << mut_ptr_copy->text << std::endl;
 *
 * // Output:
 * // Greetings from object number two
 * // Greetings from object number two
 *
 * \endcode
 */
template <typename T>
class MutablePtr : public std::shared_ptr<std::shared_ptr<T>> {
 public:
  MutablePtr(std::shared_ptr<T> ptr)
      : std::shared_ptr<std::shared_ptr<T>>(new std::shared_ptr<T>(ptr)) {}
  T& operator*() const { return **(this->get()); }
  std::shared_ptr<T> operator->() const { return *(this->get()); }
  void set(std::shared_ptr<T> ptr) { *(this->get()) = ptr; }
  explicit operator bool() const noexcept { return **this; }
};
}  // namespace common
}  // namespace everybeam
#endif
