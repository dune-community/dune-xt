// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <cassert>
#include <algorithm>

#include "pattern.hh"

namespace Dune {
namespace Stuff {
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
} // namespace Stuff
} // namespace Dune
