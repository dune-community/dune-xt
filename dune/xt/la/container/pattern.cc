// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2015 - 2016)

#include "config.h"
#include "pattern.hh"

#include <cassert>
#include <algorithm>

namespace Dune {
namespace XT {
namespace LA {

// ================================
// ==== SparsityPatternDefault ====
// ================================
SparsityPatternDefault::SparsityPatternDefault(const size_t _size)
  : vector_of_vectors_(_size)
{
}

size_t SparsityPatternDefault::size() const
{
  return vector_of_vectors_.size();
}

typename SparsityPatternDefault::InnerType& SparsityPatternDefault::inner(const size_t ii)
{
  assert(ii < size() && "Wrong index requested!");
  return vector_of_vectors_[ii];
}

const typename SparsityPatternDefault::InnerType& SparsityPatternDefault::inner(const size_t ii) const
{
  assert(ii < size() && "Wrong index requested!");
  return vector_of_vectors_[ii];
}

typename SparsityPatternDefault::ConstOuterIteratorType SparsityPatternDefault::begin() const
{
  return vector_of_vectors_.begin();
}

typename SparsityPatternDefault::ConstOuterIteratorType SparsityPatternDefault::end() const
{
  return vector_of_vectors_.end();
}

bool SparsityPatternDefault::operator==(const SparsityPatternDefault& other) const
{
  return vector_of_vectors_ == other.vector_of_vectors_;
}

bool SparsityPatternDefault::operator!=(const SparsityPatternDefault& other) const
{
  return vector_of_vectors_ != other.vector_of_vectors_;
}

SparsityPatternDefault operator+(const SparsityPatternDefault& other) const
{
  assert(other.size() == size() && "SparsityPatterns must have the same number of rows for addition!");
  SparsityPatternDefault ret = *this;
  for (size_t rr = 0; rr < size(); ++rr)
    for (const auto& cc : vector_of_vectors_[rr])
      ret.insert(rr, cc);
  ret.sort();
  return ret;
}

void SparsityPatternDefault::insert(const size_t outer_index, const size_t inner_index)
{
  assert(outer_index < size() && "Wrong index requested!");
  if (std::find(vector_of_vectors_[outer_index].begin(), vector_of_vectors_[outer_index].end(), inner_index)
      == vector_of_vectors_[outer_index].end())
    vector_of_vectors_[outer_index].push_back(inner_index);
} // ... insert(...)

void SparsityPatternDefault::sort(const size_t outer_index)
{
  assert(outer_index < size() && "Wrong index requested!");
  std::sort(vector_of_vectors_[outer_index].begin(), vector_of_vectors_[outer_index].end());
}

void SparsityPatternDefault::sort()
{
  for (auto& inner_vector : vector_of_vectors_)
    std::sort(inner_vector.begin(), inner_vector.end());
}

} // namespace LA
} // namespace XT
} // namespace Dune
