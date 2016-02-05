// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH
#define DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class IntersectionImp>
class AllNeumannBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".allneumann";
  }

  virtual const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary())
      return neumann_boundary_;
    return no_boundary_;
  }

protected:
  static constexpr NoBoundary no_boundary_{};
  static constexpr NeumannBoundary neumann_boundary_{};
}; // class AllNeumannBoundaryInfo

template <class I>
constexpr NoBoundary AllNeumannBoundaryInfo<I>::no_boundary_;
template <class I>
constexpr NeumannBoundary AllNeumannBoundaryInfo<I>::neumann_boundary_;


template <class I>
std::unique_ptr<AllNeumannBoundaryInfo<I>> make_allneumann_boundaryInfo(const Common::Configuration& /*config*/)
{
  return XT::Common::make_unique<AllNeumannBoundaryInfo<I>>();
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH
