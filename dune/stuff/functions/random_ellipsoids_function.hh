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
  DSC::FieldVector<CoordType, dim> center;
  DSC::FieldVector<CoordType, dim> radii;

  bool contains(DSC::FieldVector<CoordType, dim> point) const
  {
    const auto shifted = point - center;
    double sum = 0;
    for (auto ii : DSC::valueRange(dim)) {
      const auto pii = shifted[ii];
      sum += std::pow(pii, 2) / std::pow(radii[ii], 2);
    }
    return DSC::FloatCmp::le(sum, 1.);
  }
};

template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
          size_t rangeDimCols = 1>
class RandomEllipsoidsFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef RandomEllipsoidsFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ThisType;

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

    Localfunction(const EntityType& ent, const RangeType value)
      : BaseType(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return 0;
    }

    virtual void evaluate(const DomainType& UNUSED_UNLESS_DEBUG(xx), RangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& UNUSED_UNLESS_DEBUG(xx), JacobianRangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const RangeType value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = BaseType::dimDomain;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".RandomEllipsoidsFunction";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["lower_left"]   = "[0.0 0.0 0.0]";
    config["upper_right"]  = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"]       = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
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
    // calculate number of values and get values
    auto num_elements =
        cfg.get("num_elements", default_cfg.get<Common::FieldVector<size_t, dimDomain>>("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<RangeType> values(num_values);
    auto values_rf = cfg.get("values", default_cfg.get<std::vector<RangeFieldType>>("values"), num_values);
    for (size_t ii = 0; ii < values_rf.size(); ++ii)
      values[ii] = RangeType(values_rf[ii]);
    // create
    return Common::make_unique<ThisType>(
        cfg.get(
            "lower_left", default_cfg.get<Common::FieldVector<DomainFieldType, dimDomain>>("lower_left"), dimDomain),
        cfg.get(
            "upper_right", default_cfg.get<Common::FieldVector<DomainFieldType, dimDomain>>("upper_right"), dimDomain),
        std::move(num_elements),
        std::move(values),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  RandomEllipsoidsFunction(const Common::FieldVector<DomainFieldType, dimDomain>& lowerLeft,
                           const Common::FieldVector<DomainFieldType, dimDomain>& upperRight,
                           const Stuff::Common::Configuration& ellipsoid_cfg, const std::string nm = static_id())
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
    , name_(nm)
    , ellipsoid_cfg_(ellipsoid_cfg)
  {
    //    randomblock_config["ellipsoids.count"] = "200";
    //    randomblock_config["ellipsoids.min_radius"] = "0.01";
    //    randomblock_config["ellipsoids.max_radius"] = "0.03";
    //    randomblock_config["ellipsoids.seed"] = "0";
    //    randomblock_config["ellipsoids.children"] = "3";
    //    randomblock_config["ellipsoids.recursion_depth"] = "4";
    //    randomblock_config["ellipsoids.recursion_scale"] = "0.5";
    //    config.add(randomblock_config, "diffusion_factor");
    typedef unsigned long CT;
    const CT level_0_count = ellipsoid_cfg.get("ellipsoids.count", 200);
    const CT max_depth     = ellipsoid_cfg.get("ellipsoids.recursion_depth", 1);
    const CT children      = ellipsoid_cfg.get<CT>("ellipsoids.children", 3u, DSC::ValidateGreater<CT>(0));
    const CT total_count = level_0_count * std::pow(children, max_depth);
    ellipsoids_.resize(total_count);
    const auto seed = DomainFieldType(0);
    DSC::DefaultRNG<DomainFieldType> center_rng(
        *std::max(lowerLeft.begin(), lowerLeft.end()), *std::min(upperRight.begin(), upperRight.end()), seed);
    DSC::DefaultRNG<DomainFieldType> radii_rng(
        ellipsoid_cfg.get("ellipsoids.min_radii", 0.01), ellipsoid_cfg.get("ellipsoids.max_radii", 0.04), seed);
    for (auto ii : DSC::valueRange(0, total_count, total_count / level_0_count)) {
      std::fill(ellipsoids_[ii].center.begin(), ellipsoids_[ii].center.end(), center_rng);
      std::fill(ellipsoids_[ii].radii.begin(), ellipsoids_[ii].radii.end(), radii_rng);
    }
    for (auto current_depth : DSC::valueRange(1, max_depth + 1)) {
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
    //    std::vector< size_t > whichPartition(dimDomain, 0);
    //    const auto& ll = lowerLeft_;
    //    const auto& ur = upperRight_;
    //    const auto& ne = numElements_;
    //    for (size_t dd = 0; dd < dimDomain; ++dd) {
    //      // for points that are on upperRight_[d], this selects one partition too much
    //      // so we need to cap this
    //      whichPartition[dd] = std::min(size_t(std::floor(ne[dd]*((center[dd] - ll[dd])/(ur[dd] - ll[dd])))),
    //                                    ne[dd] - 1);
    //    }
    //    size_t subdomain = 0;
    //    if (dimDomain == 1)
    //      subdomain = whichPartition[0];
    //    else if (dimDomain == 2)
    //      subdomain = whichPartition[0] + whichPartition[1]*ne[0];
    //    else
    //      subdomain = whichPartition[0] + whichPartition[1]*ne[0] + whichPartition[2]*ne[1]*ne[0];
    //    // return the component that belongs to the subdomain
    return std::unique_ptr<Localfunction>(new Localfunction(entity, local_value));
  } // ... local_function(...)

private:
  const Common::FieldVector<DomainFieldType, dimDomain> lowerLeft_;
  const Common::FieldVector<DomainFieldType, dimDomain> upperRight_;
  std::string name_;
  const Stuff::Common::Configuration ellipsoid_cfg_;
  std::vector<Ellipsoid> ellipsoids_;
}; // class RandomEllipsoidsFunction


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_RandomEllipsoidsFunction_HH
