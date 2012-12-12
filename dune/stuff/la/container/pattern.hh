#ifndef DUNE_STUFF_LA_CONTAINER_PATTERN_HH
#define DUNE_STUFF_LA_CONTAINER_PATTERN_HH

#include <vector>
#include <set>

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {
namespace Pattern {

class Default : public std::vector<std::set<unsigned int>>
{
public:
  typedef std::vector<std::set<unsigned int>> BaseType;

  typedef BaseType::size_type size_type;

  typedef BaseType::value_type ColumnsType;

  Default(const size_type _rows)
    : BaseType(_rows)
  {
  }

  size_type rows() const
  {
    return BaseType::size();
  }

  ColumnsType& columns(const size_type _row)
  {
    assert(_row < rows() && "Wrong row requested!");
    return BaseType::operator[](_row);
  }

  const ColumnsType& columns(const size_type _row) const
  {
    assert(_row < rows() && "Wrong row requested!");
    return BaseType::operator[](_row);
  }
}; // class Default

} // namespace Pattern
} // namespace Container
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_PATTERN_HH
