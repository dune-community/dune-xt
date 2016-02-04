// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH
#define DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH

#include <string>

#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


class BoundaryType
{
protected:
  virtual std::string id() const = 0;

public:
  virtual bool operator==(const BoundaryType& other) const
  {
    return id() == other.id();
  }

  virtual bool operator!=(const BoundaryType& other) const
  {
    return !operator==(other);
  }
}; // class BoundaryType


template <class IntersectionImp>
class BoundaryInfo
{
  static_assert(is_intersection<IntersectionImp>::value, "");

public:
  typedef IntersectionImp IntersectionType;

  virtual BoundaryType type(const IntersectionType& intersection) const = 0;

  static std::string static_id()
  {
    return "xt.grid.boundaryinfo";
  }
}; // class BoundaryInfo


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH
