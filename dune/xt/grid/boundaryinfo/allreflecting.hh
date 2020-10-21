// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2016 - 2019)
//   Tobias Leibner  (2016 - 2017, 2019 - 2020)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_ALLREFLECTING_HH
#define DUNE_XT_GRID_BOUNDARYINFO_ALLREFLECTING_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include <dune/xt/grid/boundaryinfo/types.hh>

namespace Dune::XT::Grid {


static inline Common::Configuration allreflecting_boundaryinfo_default_config()
{
  return Common::Configuration({"type"}, {"xt.grid.boundaryinfo.allreflecting"});
}

// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
#endif
template <class IntersectionImp>
class AllReflectingBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  using BaseType = BoundaryInfo<IntersectionImp>;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return allreflecting_boundaryinfo_default_config().template get<std::string>("type");
  }

  const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary())
      return reflecting_boundary;
    return no_boundary;
  }
}; // class AllReflectingBoundaryInfo
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic pop
#endif

template <class I>
std::unique_ptr<AllReflectingBoundaryInfo<I>>
make_allreflecting_boundaryinfo(const Common::Configuration& /*cfg*/ = Common::Configuration())
{
  return std::make_unique<AllReflectingBoundaryInfo<I>>();
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_BOUNDARYINFO_ALLREFLECTING_HH
