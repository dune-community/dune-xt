// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_RandomEllipsoidsFunction_HH
#define DUNE_STUFF_FUNCTION_RandomEllipsoidsFunction_HH

#include <vector>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/random.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {

template <size_t dim, class CoordType = double>
struct Ellipsoid
{
  typedef DSC::FieldVector<CoordType, dim> ctype;
  ctype center;
  ctype radii;

  bool contains(ctype point) const
  {
    const auto shifted = point - center;
    double sum = 0;
    for (auto ii : DSC::valueRange(dim)) {
      const auto pii = shifted[ii];
      sum += std::pow(pii, 2) / std::pow(radii[ii], 2);
    }
    return DSC::FloatCmp::le(sum, 1.);
  }
  bool intersects_cube(ctype ll, ctype ur) const
  {
    DUNE_THROW(NotImplemented, "");
  }
};

template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
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
      , value_(value)
      , local_ellipsoids_(local_ellipsoids)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      for (const auto& ellipsoid : local_ellipsoids_) {
        if (ellipsoid.contains(xx)) {
          ret = value_;
          return;
        }
      }
      ret = 0;
    }

    virtual void jacobian(const DomainType& UNUSED_UNLESS_DEBUG(xx), JacobianRangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
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
    config["lower_left"]  = "[0.0 0.0 0.0]";
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
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    DUNE_THROW(NotImplemented, "");
  } // ... create(...)

  RandomEllipsoidsFunction(const Common::FieldVector<DomainFieldType, dimDomain>& lowerLeft,
                           const Common::FieldVector<DomainFieldType, dimDomain>& upperRight,
                           const Stuff::Common::Configuration& ellipsoid_cfg, const std::string nm = static_id())
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
    , name_(nm)
    , ellipsoid_cfg_(ellipsoid_cfg)
  {
    //    randomblock_config["ellipsoids.seed"] = "0";
    //    randomblock_config["ellipsoids.children"] = "3";
    //    randomblock_config["ellipsoids.recursion_depth"] = "4";
    //    randomblock_config["ellipsoids.recursion_scale"] = "0.5";
    //    config.add(randomblock_config, "diffusion_factor");
    typedef unsigned long CT;
    const CT level_0_count = ellipsoid_cfg.get("ellipsoids.count", 10);
    const CT max_depth     = ellipsoid_cfg.get("ellipsoids.recursion_depth", 1);
    const CT children      = ellipsoid_cfg.get<CT>("ellipsoids.children", 3u, DSC::ValidateGreater<CT>(0));
    const CT total_count = level_0_count * std::pow(children, max_depth);
    ellipsoids_.resize(level_0_count);
    ellipsoids_.reserve(total_count);
    const auto seed = DomainFieldType(ellipsoid_cfg.get("ellipsoids.seed", 0));
    DSC::DefaultRNG<DomainFieldType> center_rng(
        *std::max(lowerLeft.begin(), lowerLeft.end()), *std::min(upperRight.begin(), upperRight.end()), seed);
    DSC::DefaultRNG<DomainFieldType> radii_rng(
        ellipsoid_cfg.get("ellipsoids.min_radius", 0.01), ellipsoid_cfg.get("ellipsoids.max_radius", 0.04), seed);
    DSC::DefaultRNG<DomainFieldType> dist_rng(-ellipsoid_cfg.get("ellipsoids.max_radius", 0.04) / 2.,
                                              ellipsoid_cfg.get("ellipsoids.max_radius", 0.04) / 2.,
                                              seed);

    const auto parent_range = DSC::valueRange(0ul, level_0_count);
    for (auto ii : parent_range) {
      std::generate(ellipsoids_[ii].center.begin(), ellipsoids_[ii].center.end(), center_rng);
      std::generate(ellipsoids_[ii].radii.begin(), ellipsoids_[ii].radii.end(), radii_rng);
    }

    std::function<void(CT, const EllipsoidType&)> recurse_add = [&](CT current_level, const EllipsoidType& parent) {
      if (current_level > max_depth)
        return;
      for (auto child_no : DSC::valueRange(children)) {
        EllipsoidType child = parent;
        auto displacement = [&](DomainFieldType& coord) {
          coord += dist_rng() * std::pow(ellipsoid_cfg.get("ellipsoids.recursion_scale", 0.5), current_level);
        };
        std::for_each(child.center.begin(), child.center.end(), displacement);
        std::generate(child.radii.begin(), child.radii.end(), radii_rng);
        ellipsoids_.push_back(child);
        recurse_add(current_level + 1, child);
      }
    };

    for (auto ii : parent_range) {
      recurse_add(0, ellipsoids_[ii]);
    }
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

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    constexpr auto local_value = DomainFieldType(1);
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    bool in_ellipsoid = false;
    std::vector<EllipsoidType> local_ellipsoids;
    for (const auto& ellipsoid : ellipsoids_) {
      typename EllipsoidType::ctype ll, ur;
      if (ellipsoid.intersects_cube(ll, ur)) {
        local_ellipsoids.push_back(ellipsoid);
      }
    }
    return std::unique_ptr<Localfunction>(new Localfunction(entity, local_value, std::move(local_ellipsoids)));
  } // ... local_function(...)

private:
  const Common::FieldVector<DomainFieldType, dimDomain> lowerLeft_;
  const Common::FieldVector<DomainFieldType, dimDomain> upperRight_;
  std::string name_;
  const Stuff::Common::Configuration ellipsoid_cfg_;
  std::vector<EllipsoidType> ellipsoids_;
}; // class RandomEllipsoidsFunction


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_RandomEllipsoidsFunction_HH
