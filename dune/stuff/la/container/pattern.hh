// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_CONTAINER_PATTERN_HH
#define DUNE_STUFF_LA_CONTAINER_PATTERN_HH

#include <cstddef>
#include <vector>
#include <set>

namespace Dune {
namespace Stuff {
namespace LA {

class SparsityPatternDefault
{
private:
  typedef std::vector<std::vector<size_t>> BaseType;

public:
  typedef BaseType::value_type InnerType;
  typedef typename BaseType::const_iterator ConstOuterIteratorType;

  explicit SparsityPatternDefault(const size_t _size = 0);

  size_t size() const;

  InnerType& inner(const size_t ii);

  const InnerType& inner(const size_t ii) const;

  ConstOuterIteratorType begin() const;

  ConstOuterIteratorType end() const;

  bool operator==(const SparsityPatternDefault& other) const;

  bool operator!=(const SparsityPatternDefault& other) const;

  void insert(const size_t outer_index, const size_t inner_index);

  void sort(const size_t outer_index);

  void sort();

private:
  BaseType vector_of_vectors_;
}; // class SparsityPatternDefault

} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_PATTERN_HH
