// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
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
template <class IntersectionImp>
class NormalBasedBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;
  typedef NormalBasedBoundaryInfo<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::WorldType;
  using BaseType::dimDomain;
  using BaseType::dimWorld;

  static std::string static_id()
  {
    return normalbased_boundaryinfo_default_config().template get<std::string>("type");
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration cfg = normalbased_boundaryinfo_default_config())
  {
    const Common::Configuration default_cfg = normalbased_boundaryinfo_default_config();
    // get default
    const std::string default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    if (default_type != "dirichlet" && default_type != "neumann")
      DUNE_THROW(Common::Exceptions::configuration_error, "Wrong default '" << default_type << "' given!");
    const bool default_to_dirichlet = default_type == "dirichlet";
    // get tolerance
    const DomainFieldType tol = cfg.get("compare_tolerance", default_cfg.get<DomainFieldType>("compare_tolerance"));
    // get dirichlet and neumann
    std::vector<WorldType> dirichlets = getVectors(cfg, "dirichlet");
    std::vector<WorldType> neumanns = getVectors(cfg, "neumann");
    // return
    return Common::make_unique<ThisType>(default_to_dirichlet, dirichlets, neumanns, tol);
  } // ... create(...)

  NormalBasedBoundaryInfo(const bool default_to_dirichlet = true,
                          const std::vector<WorldType> dirichlet_normals = std::vector<WorldType>(),
                          const std::vector<WorldType> neumann_normals = std::vector<WorldType>(),
                          const DomainFieldType tol = 1e-10)
    : default_to_dirichlet_(default_to_dirichlet)
    , dirichlet_normals_(dirichlet_normals)
    , neumann_normals_(neumann_normals)
    , tol_(tol)
  {
    // normalize
    for (auto& normal : dirichlet_normals_)
      normal /= normal.two_norm();
    for (auto& normal : neumann_normals_)
      normal /= normal.two_norm();
    // sanity check
    for (auto& dirichletNormal : dirichlet_normals_) {
      if (contains(dirichletNormal, neumann_normals_))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "Given normals are too close for given tolerance '" << tol << "'!");
    }
  } // NormalBased(...)

  NormalBasedBoundaryInfo(const std::vector<WorldType> dirichlet_normals,
                          const std::vector<WorldType> neumann_normals,
                          const std::vector<WorldType> reflecting_normals,
                          const DomainFieldType tol = 1e-10)
    : NormalBasedBoundaryInfo(true, dirichlet_normals, neumann_normals, tol)
  {
    reflecting_normals_ = reflecting_normals;
    // normalize
    for (auto& normal : reflecting_normals_)
      normal /= normal.two_norm();
    // sanity checks
    for (auto& dirichletNormal : dirichlet_normals_) {
      if (contains(dirichletNormal, reflecting_normals_))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "Given normals are too close for given tolerance '" << tol << "'!");
    }
    for (auto& neumannNormal : neumann_normals_) {
      if (contains(neumannNormal, reflecting_normals_))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "Given normals are too close for given tolerance '" << tol << "'!");
    }
  } // NormalBased(...)


  virtual const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return no_boundary_;
    const WorldType outerNormal = intersection.centerUnitOuterNormal();
    if (contains(outerNormal, dirichlet_normals_))
      return dirichlet_boundary_;
    else if (contains(outerNormal, neumann_normals_))
      return neumann_boundary_;
    else if (contains(outerNormal, reflecting_normals_))
      return reflecting_boundary_;
    else
      return dirichlet_boundary_;
  } // ... type(...)

protected:
  static std::vector<WorldType> getVectors(const Common::Configuration& config, const std::string key)
  {
    std::vector<WorldType> ret;
    if (config.has_sub(key)) {
      bool found = true;
      size_t counter = 0;
      while (found) {
        const std::string localKey = key + "." + Dune::XT::Common::to_string(counter);
        if (config.has_key(localKey))
          ret.push_back(config.get<WorldType>(localKey, dimWorld));
        else
          found = false;
        ++counter;
      }
    } else if (config.has_key(key))
      ret.push_back(config.get<WorldType>(key, dimWorld));
    return ret;
  } // ... getVectors(...)

  bool contains(const WorldType& normal, const std::vector<WorldType>& vectors) const
  {
    for (auto& vector : vectors)
      if (Dune::XT::Common::FloatCmp::eq(normal, vector, tol_))
        return true;
    return false;
  }

  const bool default_to_dirichlet_;
  std::vector<WorldType> dirichlet_normals_;
  std::vector<WorldType> neumann_normals_;
  std::vector<WorldType> reflecting_normals_;
  const DomainFieldType tol_;
  static constexpr NoBoundary no_boundary_{};
  static constexpr DirichletBoundary dirichlet_boundary_{};
  static constexpr NeumannBoundary neumann_boundary_{};
  static constexpr ReflectingBoundary reflecting_boundary_{};
}; // class NormalBasedBoundaryInfo
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic pop
#endif

template <class I>
constexpr NoBoundary NormalBasedBoundaryInfo<I>::no_boundary_;
template <class I>
constexpr DirichletBoundary NormalBasedBoundaryInfo<I>::dirichlet_boundary_;
template <class I>
constexpr NeumannBoundary NormalBasedBoundaryInfo<I>::neumann_boundary_;
template <class I>
constexpr ReflectingBoundary NormalBasedBoundaryInfo<I>::reflecting_boundary_;


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
