// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014 - 2015)

#ifndef DUNE_XT_FUNCTIONS_RANDOMELLIPSOIDS_HH
#define DUNE_XT_FUNCTIONS_RANDOMELLIPSOIDS_HH

#include <cmath>
#include <memory>
#include <vector>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/debug.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/random.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {

template <size_t dim, class CoordType = double>
struct Ellipsoid
{
  typedef FieldVector<CoordType, dim> DomainType;
  DomainType center;
  DomainType radii;

  bool contains(DomainType point) const
  {
    const auto shifted = point - center;
    double sum = 0;
    for (auto ii : Common::value_range(dim)) {
      const auto pii = shifted[ii];
      sum += std::pow(pii, 2) / std::pow(radii[ii], 2);
    }
    return FloatCmp::le(sum, 1.);
  }
  bool intersects_cube(DomainType /*ll*/, DomainType /*ur*/) const
  {
    DUNE_THROW(NotImplemented, "");
  }
};

template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class RandomEllipsoidsFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
protected:
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef RandomEllipsoidsFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ThisType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;
  typedef Ellipsoid<dimDomain, DomainFieldType> EllipsoidType;
  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeFieldType RangeFieldType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const RangeType value, std::vector<EllipsoidType>&& local_ellipsoids)
      : BaseType(ent)
      , geometry_(ent.geometry())
      , value_(value)
      , local_ellipsoids_(local_ellipsoids)
    {
      //      DXTC_LOG_DEBUG_0 << "create local LF Ellips with " << local_ellipsoids_.size() << " instances\n";
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
    {
      return 8;
    }

    virtual void evaluate(const DomainType& xx_local,
                          RangeType& ret,
                          const Common::Parameter& /*mu*/ = Common::Parameter()) const override
    {
      assert(this->is_a_valid_point(xx_local));
      const auto xx_global = geometry_.global(xx_local);
      for (const auto& ellipsoid : local_ellipsoids_) {
        if (ellipsoid.contains(xx_global)) {
          ret = value_;
          //          DXTC_LOG_DEBUG_0 << "ell  INSIDE " << ellipsoid.center << " with xx " << xx_global << "\n";
          return;
        }
      }
      ret = 0;
    }

    virtual void jacobian(const DomainType& DXTC_DEBUG_ONLY(xx),
                          JacobianRangeType& ret,
                          const Common::Parameter& /*mu*/ = Common::Parameter()) const override
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const typename EntityImp::Geometry geometry_;
    const RangeType value_;
    const std::vector<EllipsoidType> local_ellipsoids_;
  }; // class Localfunction

public:
  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".RandomEllipsoidsFunction";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    DUNE_THROW(NotImplemented, "");
  } // ... create(...)

  RandomEllipsoidsFunction(const Common::FieldVector<DomainFieldType, dimDomain>& lowerLeft,
                           const Common::FieldVector<DomainFieldType, dimDomain>& upperRight,
                           const Common::Configuration& ellipsoid_cfg,
                           const std::string nm = static_id())
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
    , name_(nm)
    , ellipsoid_cfg_(ellipsoid_cfg)
  {
    typedef unsigned long UL;
    const UL level_0_count = ellipsoid_cfg.get("ellipsoids.count", 10);
    const UL max_depth = ellipsoid_cfg.get("ellipsoids.recursion_depth", 1);
    const UL children = ellipsoid_cfg.get<UL>("ellipsoids.children", 3u); //, ValidateLess<UL>(0));
    const UL total_count = level_0_count + level_0_count * std::pow(children, max_depth + 1);
    ellipsoids_.resize(level_0_count);
    ellipsoids_.reserve(total_count);
    const auto seed = DomainFieldType(ellipsoid_cfg.get("ellipsoids.seed", 0));
    const auto min_radius = ellipsoid_cfg.get("ellipsoids.min_radius", 0.01);
    const auto max_radius = ellipsoid_cfg.get("ellipsoids.max_radius", 0.02);
    const auto child_displacement = ellipsoid_cfg.get("ellipsoids.max_child_displacement", max_radius);
    typedef Common::DefaultRNG<DomainFieldType> RNG;
    RNG center_rng(0, 1, seed);
    RNG radii_rng(min_radius, max_radius, seed);
    RNG dist_rng(min_radius, child_displacement, seed);
    RNG sign_rng(-1, 1, seed);

    const auto parent_range = Common::value_range(0ul, level_0_count);
    for (auto ii : parent_range) {
      std::generate(
          ellipsoids_[ii].center.begin(), ellipsoids_[ii].center.end(), [&center_rng]() { return center_rng(); });
      std::generate(ellipsoids_[ii].radii.begin(), ellipsoids_[ii].radii.end(), [&radii_rng]() { return radii_rng(); });
    }

    std::function<void(UL, const EllipsoidType&)> recurse_add = [&](UL current_level, const EllipsoidType& parent) {
      if (current_level > max_depth)
        return;
      for (const auto unused_counter : Common::value_range(children)) {
        EllipsoidType child = parent;
        const double scale = std::pow(ellipsoid_cfg.get("ellipsoids.recursion_scale", 0.5), current_level);
        const auto displace = [&](DomainFieldType& coord) {
          const auto disp = dist_rng() * signum(sign_rng());
          coord += disp;
          DXTC_LOG_DEBUG_0 << disp << ";";
        };
        std::for_each(child.center.begin(), child.center.end(), displace);
        std::generate(child.radii.begin(), child.radii.end(), [&radii_rng, scale]() { return radii_rng() * scale; });
        ellipsoids_.push_back(child);
        recurse_add(current_level + 1, child);
      }
    };

    for (auto ii : parent_range) {
      recurse_add(0, ellipsoids_[ii]);
    }
    DXTC_LOG_DEBUG_0 << "generated " << ellipsoids_.size() << " of " << total_count << "\n";
    to_file(*Common::make_ofstream("ellipsoids.txt"));
  }

  RandomEllipsoidsFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".RandomEllipsoidsFunction";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  void to_file(std::ostream& out) const
  {
    boost::format line("%d|%g;%g|%g;%g\n");
    size_t id = 0;
    for (const auto& ellipsoid : ellipsoids_) {
      const auto str =
          (line % id % ellipsoid.center[0] % ellipsoid.center[1] % ellipsoid.radii[0] % ellipsoid.radii[1]).str();
      out << str;
    }
    out.flush();
  }

private:
  std::tuple<typename EllipsoidType::DomainType, typename EllipsoidType::DomainType>
  bounding_box(const EntityType& entity) const
  {
    typename EllipsoidType::DomainType ll, ur;
    typedef Common::MinMaxAvg<DomainFieldType> MinMaxAvgType;
    std::array<MinMaxAvgType, dimDomain> coord_limits;
    const auto& geo = entity.geometry();
    for (auto i : Common::value_range(geo.corners())) {
      const auto& corner(geo.corner(i));
      for (size_t k = 0; k < dimDomain; ++k)
        coord_limits[k](corner[k]);
    }
    for (auto ii : Common::value_range(dimDomain)) {
      ll[ii] = coord_limits[ii].min();
      ur[ii] = coord_limits[ii].max();

      return std::make_pair(std::move(ll), std::move(ur));
    }
  }

public:
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    const auto local_value = DomainFieldType(ellipsoid_cfg_.get("ellipsoids.local_value", 1.));
    // decide on the subdomain the center of the entity belongs to
    const auto max_radius = ellipsoid_cfg_.get("ellipsoids.max_radius", 0.04);
    typename EllipsoidType::DomainType ll, ur;
    std::tie(ll, ur) = bounding_box(entity);
    std::for_each(ll.begin(), ll.end(), [&max_radius](DomainFieldType& lc) { lc -= 1.1 * max_radius; });
    std::for_each(ur.begin(), ur.end(), [&max_radius](DomainFieldType& uc) { uc += 1.1 * max_radius; });
    std::vector<EllipsoidType> local_ellipsoids;
    for (const auto& ellipsoid : ellipsoids_) {
      bool inside = true;
      for (auto ii : Common::value_range(dimDomain)) {
        inside = inside && FloatCmp::ge(ellipsoid.center[ii], ll[ii]);
        inside = inside && FloatCmp::le(ellipsoid.center[ii], ur[ii]);
      }
      if (inside) {
        local_ellipsoids.push_back(ellipsoid);
      }
    }
    std::vector<EllipsoidType> tmp = ellipsoids_;
    return std::unique_ptr<Localfunction>(new Localfunction(entity, local_value, std::move(tmp)));
  } // ... local_function(...)

private:
  const Common::FieldVector<DomainFieldType, dimDomain> lowerLeft_;
  const Common::FieldVector<DomainFieldType, dimDomain> upperRight_;
  const std::string name_;
  const Common::Configuration ellipsoid_cfg_;
  std::vector<EllipsoidType> ellipsoids_;
}; // class RandomEllipsoidsFunction

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_RANDOMELLIPSOIDS_HH
