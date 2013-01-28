#ifndef DUNE_STUFF_LA_CONTAINER_PATTERN_HH
#define DUNE_STUFF_LA_CONTAINER_PATTERN_HH

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

  typedef BaseType::value_type SetType;

  SparsityPatternDefault(const size_type _size)
    : vectorOfSets_(_size)
  {
  }

  size_type size() const
  {
    return vectorOfSets_.size();
  }

  SetType& set(const size_type _index)
  {
    assert(_index < size() && "Wrong index requested!");
    return vectorOfSets_[_index];
  }

  const SetType& set(const size_type _index) const
  {
    assert(_index < size() && "Wrong index requested!");
    return vectorOfSets_[_index];
  }

private:
  BaseType vectorOfSets_;
}; // class SparsityPatternDefault

} // namespace Container
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_PATTERN_HH
