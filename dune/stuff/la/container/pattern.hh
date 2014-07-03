// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
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

  SparsityPatternDefault(const size_t _size = 0);

  size_t size() const;

  InnerType& inner(const size_t ii);

  const InnerType& inner(const size_t ii) const;

  ConstOuterIteratorType begin() const;

  ConstOuterIteratorType end() const;

  bool operator==(const SparsityPatternDefault& other) const;

  bool operator!=(const SparsityPatternDefault& other) const;

private:
  BaseType vectorOfSets_;
}; // class SparsityPatternDefault


//// forward of the Interface
// template< class T >
// class MatrixInterface;


// template< class T >
// SparsityPatternDefault* createCompressedSparsityPattern(const SparsityPatternDefault& uncompressedPattern,
//                                                        const MatrixInterface< T >& matrix,
//                                                        const typename T::ElementType threshhold = typename
//                                                        T::ElementType(0))
//{
//  assert(!(threshhold < typename T::ElementType(0)) && "Please provide a nonnegative threshhold!");
//  // we create a new pattern of only nonzero entries
//  SparsityPatternDefault* compressedPattern = new SparsityPatternDefault(uncompressedPattern.size());
//  typename T::ElementType absValue(0);
//  // * therefore we loop over the uncompressed pattern
//  for (typename SparsityPatternDefault::size_t row = 0; row < uncompressedPattern.size(); ++row) {
//    // * get the uncompressed row,
//    const auto& uncompressedRowSet = uncompressedPattern.set(row);
//    // * get the new one
//    auto& compressedRowSet = compressedPattern->set(row);
//    // * and loop over the uncompressed row
//    for (auto col : uncompressedRowSet) {
//      // * get the value in the matric
//      absValue = std::abs(matrix.get(row, col));
//      // * and add the column to the new pattern, if the value is large enough
//      if (absValue > threshhold)
//        compressedRowSet.insert(col);
//    }
//  }
//  return compressedPattern;
//} // ... compressSparsityPattern(...)


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_PATTERN_HH
