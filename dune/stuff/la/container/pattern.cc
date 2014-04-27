// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "pattern.hh"

#include <assert.h>

namespace Dune {
namespace Stuff {
namespace LA {


// ================================
// ==== SparsityPatternDefault ====
// ================================
SparsityPatternDefault::SparsityPatternDefault(const size_t _size)
  : vectorOfSets_(_size)
{
}

size_t SparsityPatternDefault::size() const
{
  return vectorOfSets_.size();
}

typename SparsityPatternDefault::InnerType& SparsityPatternDefault::inner(const size_t ii)
{
  assert(ii < size() && "Wrong index requested!");
  return vectorOfSets_[ii];
}

const typename SparsityPatternDefault::InnerType& SparsityPatternDefault::inner(const size_t ii) const
{
  assert(ii < size() && "Wrong index requested!");
  return vectorOfSets_[ii];
}

typename SparsityPatternDefault::ConstOuterIteratorType SparsityPatternDefault::begin() const
{
  return vectorOfSets_.begin();
}

typename SparsityPatternDefault::ConstOuterIteratorType SparsityPatternDefault::end() const
{
  return vectorOfSets_.end();
}

bool SparsityPatternDefault::operator==(const SparsityPatternDefault& other) const
{
  return vectorOfSets_ == other.vectorOfSets_;
}

bool SparsityPatternDefault::operator!=(const SparsityPatternDefault& other) const
{
  return vectorOfSets_ != other.vectorOfSets_;
}


} // namespace LA
} // namespace Stuff
} // namespace Dune
