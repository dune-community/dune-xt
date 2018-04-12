// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH
#define DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH

#include <unordered_map>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/exceptions.hh>

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


static inline Common::Configuration normalbased_boundaryinfo_default_config()
{
  Common::Configuration config;
  config["type"] = "xt.grid.boundaryinfo.normalbased";
  config["default"] = NoBoundary().id();
  config["compare_tolerance"] = "1e-10";
  config["0." + DirichletBoundary().id()] = "[1 0 0 0]";
  config["1." + NeumannBoundary().id()] = "[0 1 0 0]";
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

  static std::unique_ptr<ThisType> create(const Common::Configuration& cfg = normalbased_boundaryinfo_default_config())
  {
    const Common::Configuration default_cfg = normalbased_boundaryinfo_default_config();
    // get tolerance and default
    const auto tol = cfg.get("compare_tolerance", default_cfg.get<DomainFieldType>("compare_tolerance"));
    const auto default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    // create
    auto ret = std::make_unique<ThisType>(tol, make_boundary_type(default_type));
    // get other normals and boundary types and register
    bool search_for_types = true;
    size_t counter = 0;
    std::unique_ptr<BoundaryType> boundary_type(nullptr);
    WorldType normal;
    while (search_for_types) {
      if (cfg.has_sub(Common::to_string(counter))) {
        const Common::Configuration sub_cfg = cfg.sub(Common::to_string(counter));
        if (sub_cfg.getValueKeys().size() != 1)
          DUNE_THROW(Exceptions::boundary_info_error,
                     "while processing sub config "
                         << counter
                         << " of cfg (see below): could not parse give config."
                         << "For each normal, you have to provide a sub config with exactly "
                            "one key/value pair, where the key determines the BoundaryType "
                            "and the value determines the normal (see below for a valid default config)."
                         << "\n\n   This was the given config:\n"
                         << cfg
                         << "\n\n   This is a suitable default config:\n"
                         << default_cfg);
        const auto boundary_type_key = sub_cfg.getValueKeys()[0];
        try {
          boundary_type = std::unique_ptr<BoundaryType>(make_boundary_type(boundary_type_key));
        } catch (const Exceptions::boundary_type_error& ee) {
          DUNE_THROW(
              Exceptions::boundary_info_error,
              "while processing sub config "
                  << counter
                  << " of cfg (see below): given key is not a valid BoundaryType (see below for the original error)."
                  << "For each normal, you have to provide a sub config with exactly "
                     "one key/value pair, where the key determines the BoundaryType "
                     "and the value determines the normal (see below for a valid default config)."
                  << "\n\n   This was the given config:\n"
                  << cfg
                  << "\n\n   This is a suitable default config:\n"
                  << default_cfg
                  << "\nn   This was the original error:\n"
                  << ee.what());
        }
        try {
          normal = sub_cfg.get<WorldType>(boundary_type_key);
        } catch (const Common::Exceptions::configuration_error& ee) {
          DUNE_THROW(Exceptions::boundary_info_error,
                     "while processing sub config "
                         << counter
                         << " of cfg (see below): given value is not a valid normal (see below for the original error)."
                         << "For each normal, you have to provide a sub config with exactly "
                            "one key/value pair, where the key determines the BoundaryType "
                            "and the value determines the normal (see below for a valid default config)."
                         << "\n\n   This was the given config:\n"
                         << cfg
                         << "\n\n   This is a suitable default config:\n"
                         << default_cfg
                         << "\nn   This was the original error:\n"
                         << ee.what());
        }
        ret->register_new_normal(normal, boundary_type->copy());
        ++counter;
      } else
        search_for_types = false;
    }
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
    normal_to_type_map_.emplace(normal, std::shared_ptr<BoundaryType>(std::move(boundary_type)));
  } // ... void register_new_normal(...)

  const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return no_boundary;
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
  std::unordered_map<WorldType, std::shared_ptr<BoundaryType>> normal_to_type_map_;
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
