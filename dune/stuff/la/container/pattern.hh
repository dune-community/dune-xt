#ifndef DUNE_STUFF_LA_CONTAINER_PATTERN_HH
#define DUNE_STUFF_LA_CONTAINER_PATTERN_HH

#include <cmath>
#include <vector>
#include <set>
#include <assert.h>

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {


class SparsityPatternDefault
{
private:
  typedef std::vector<std::set<unsigned int>> BaseType;

public:
  typedef BaseType::size_type size_type;
  typedef BaseType::value_type InnerType;
  typedef typename BaseType::const_iterator ConstOuterIteratorType;

  SparsityPatternDefault(const size_type _size)
    : vectorOfSets_(_size)
  {
  }

  size_type size() const
  {
    return vectorOfSets_.size();
  }

  InnerType& inner(const size_type ii)
  {
    assert(ii < size() && "Wrong index requested!");
    return vectorOfSets_[ii];
  }

  const InnerType& inner(const size_type ii) const
  {
    assert(ii < size() && "Wrong index requested!");
    return vectorOfSets_[ii];
  }

  ConstOuterIteratorType begin() const
  {
    return vectorOfSets_.begin();
  }

  ConstOuterIteratorType end() const
  {
    return vectorOfSets_.end();
  }

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
//  for (typename SparsityPatternDefault::size_type row = 0; row < uncompressedPattern.size(); ++row) {
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


} // namespace Container
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_PATTERN_HH
