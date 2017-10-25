// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH
#define DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


static inline Common::Configuration normalbased_boundaryinfo_default_config()
{
  Common::Configuration config;
  config["type"] = "xt.grid.boundaryinfo.normalbased";
  config["default"] = "dirichlet";
  config["compare_tolerance"] = "1e-10";
  config["neumann.0"] = "[1 0 0 0]";
  return config;
}

// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif

/**
 * Use as in
\code
XT::Grid::NormalBasedBoundaryInfo<...> boundary_info;
boundary_info.register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
boundary_info.register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
\endcode
 **/
template <class I>
class NormalBasedBoundaryInfo : public BoundaryInfo<I>
{
  using BaseType = BoundaryInfo<I>;
  using ThisType = NormalBasedBoundaryInfo<I>;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::IntersectionType;
  using typename BaseType::WorldType;

  static std::string static_id()
  {
    return normalbased_boundaryinfo_default_config().template get<std::string>("type");
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration cfg = normalbased_boundaryinfo_default_config())
  {
    DUNE_THROW(NotImplemented, "Until we have a factory for BoundaryTypes!");
    const Common::Configuration default_cfg = normalbased_boundaryinfo_default_config();
    // get tolerance and default
    const DomainFieldType tol = cfg.get("compare_tolerance", default_cfg.get<DomainFieldType>("compare_tolerance"));
    // ...
    // create
    auto ret = std::make_unique<ThisType>(tol /*, default=*/);
    // get other normals and boundary types and register
    // ...
    // return
    return ret;
  } // ... create(...)

  /**
   * \attention Takes ownership of default_boundary_type, do not delete manually!
   */
  NormalBasedBoundaryInfo(const DomainFieldType tol = 1e-10, BoundaryType*&& default_boundary_type = new NoBoundary())
    : tol_(tol)
    , default_boundary_type_(std::move(default_boundary_type))
  {
  }

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  void register_new_normal(const WorldType& normal, BoundaryType*&& boundary_type)
  {
    const auto normalized_normal = normal / normal.two_norm();
    for (const auto& normal_and_type_pair : normal_to_type_map_) {
      const auto& existing_normal = normal_and_type_pair.first;
      if (XT::Common::FloatCmp::eq(existing_normal, normalized_normal, tol_))
        DUNE_THROW(InvalidStateException, "Given normals are too close for given tolerance '" << tol_ << "'!");
    }
    normal_to_type_map_.emplace(normal, std::move(boundary_type));
  } // ... void register_new_normal(...)

  const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return no_boundary_;
    const WorldType outer_normal = intersection.centerUnitOuterNormal();
    for (const auto& normal_and_type_pair : normal_to_type_map_) {
      const auto& normal = normal_and_type_pair.first;
      const auto& type_ptr = normal_and_type_pair.second;
      if (XT::Common::FloatCmp::eq(outer_normal, normal, tol_))
        return *type_ptr;
    }
    return *default_boundary_type_;
  } // ... type(...)

private:
  const DomainFieldType tol_;
  const std::unique_ptr<BoundaryType> default_boundary_type_;
  const NoBoundary no_boundary_;
  std::map<WorldType, std::shared_ptr<BoundaryType>> normal_to_type_map_;
}; // class NormalBasedBoundaryInfo
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic pop
#endif


template <class I>
std::unique_ptr<NormalBasedBoundaryInfo<I>>
make_normalbased_boundaryinfo(const Common::Configuration& cfg = normalbased_boundaryinfo_default_config())
{
  return NormalBasedBoundaryInfo<I>::create(cfg);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH
