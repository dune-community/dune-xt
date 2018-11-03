// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH
#define DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include <dune/xt/grid/boundaryinfo/types.hh>

namespace Dune {
namespace XT {
namespace Grid {


static inline Common::Configuration alldirichlet_boundaryinfo_default_config()
{
  return Common::Configuration({"type"}, {"xt.grid.boundaryinfo.alldirichlet"});
}

// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif
template <class IntersectionImp>
class AllDirichletBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return alldirichlet_boundaryinfo_default_config().template get<std::string>("type");
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
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic pop
#endif

template <class I>
constexpr NoBoundary AllDirichletBoundaryInfo<I>::no_boundary_;
template <class I>
constexpr DirichletBoundary AllDirichletBoundaryInfo<I>::dirichlet_boundary_;


template <class I>
std::unique_ptr<AllDirichletBoundaryInfo<I>>
make_alldirichlet_boundaryinfo(const Common::Configuration& /*cfg*/ = Common::Configuration())
{
  return XT::Common::make_unique<AllDirichletBoundaryInfo<I>>();
}

template <class GV>
std::unique_ptr<AllDirichletBoundaryInfo<extract_intersection_t<GV>>>
make_alldirichlet_boundaryinfo(GV, const Common::Configuration& /*cfg*/ = Common::Configuration())
{
  return XT::Common::make_unique<AllDirichletBoundaryInfo<extract_intersection_t<GV>>>();
}

template <class IntersectionImp>
class ProcessBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return "xt.grid.boundaryinfo.process";
  }

  virtual const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.neighbor() && !intersection.boundary())
      return dirichlet_boundary_;
    return no_boundary_;
  }

protected:
  static constexpr NoBoundary no_boundary_{};
  static constexpr DirichletBoundary dirichlet_boundary_{};
}; // class ProcessBoundaryInfo
template <class I>
constexpr NoBoundary ProcessBoundaryInfo<I>::no_boundary_;
template <class I>
constexpr DirichletBoundary ProcessBoundaryInfo<I>::dirichlet_boundary_;

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_ALLDIRICHLET_HH
