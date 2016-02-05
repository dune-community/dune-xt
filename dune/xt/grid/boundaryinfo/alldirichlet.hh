// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH
#define DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class IntersectionImp>
class AllDirichletBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".alldirichlet";
  }

  virtual const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary())
      return dirichlet_boundary_;
    return no_boundary_;
  }

protected:
  static constexpr NoBoundary no_boundary_{};
  static constexpr DirichletBoundary dirichlet_boundary_{};
}; // class AllDirichletBoundaryInfo

template <class I>
constexpr NoBoundary AllDirichletBoundaryInfo<I>::no_boundary_;
template <class I>
constexpr DirichletBoundary AllDirichletBoundaryInfo<I>::dirichlet_boundary_;

Common::Configuration alldirichlet_boundaryinfo_default_config()
{
  return Common::Configuration("type", "xt.grid.boundaryinfo.alldirichlet");
}


template <class I>
std::unique_ptr<AllDirichletBoundaryInfo<I>>
make_alldirichlet_boundaryinfo(const Common::Configuration& /*cfg*/ = Common::Configuration())
{
  return XT::Common::make_unique<AllDirichletBoundaryInfo<I>>();
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH
