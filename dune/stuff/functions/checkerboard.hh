// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/debug.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Checkerboard
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef Checkerboard<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange     = BaseType::dimRange;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const RangeType value)
      : BaseType(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE
    {
      return 0;
    }

    virtual void evaluate(const DomainType& UNUSED_UNLESS_DEBUG(xx), RangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& UNUSED_UNLESS_DEBUG(xx), JacobianRangeType& ret) const DS_OVERRIDE
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
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange     = BaseType::dimRange;
  static const int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".checkerboard";
  }

  static Common::ConfigContainer default_config(const std::string sub_name = "")
  {
    Common::ConfigContainer config;
    config["lower_left"]   = "[0.0 0.0 0.0]";
    config["upper_right"]  = "[1.0 1.0 1.0]";
    config["num_elements"] = "[2 2 2]";
    config["values"]       = "[1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigContainer tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigContainer config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::ConfigContainer cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::ConfigContainer default_cfg = default_config();
    // calculate number of values and get values
    auto num_elements = cfg.get("num_elements", default_cfg.get<std::vector<size_t>>("num_elements"), dimDomain);
    size_t num_values = 1;
    for (size_t ii = 0; ii < num_elements.size(); ++ii)
      num_values *= num_elements[ii];
    std::vector<RangeType> values(num_values);
    auto values_rf = cfg.get("values", default_cfg.get<std::vector<RangeFieldType>>("values"), num_values);
    for (size_t ii = 0; ii < values_rf.size(); ++ii)
      values[ii] = RangeType(values_rf[ii]);
    // create
    return Common::make_unique<ThisType>(
        cfg.get("lower_left", default_cfg.get<std::vector<DomainFieldType>>("lower_left"), dimDomain),
        cfg.get("upper_right", default_cfg.get<std::vector<DomainFieldType>>("upper_right"), dimDomain),
        std::move(num_elements),
        std::move(values),
        cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Checkerboard(std::vector<DomainFieldType>&& lowerLeft, std::vector<DomainFieldType>&& upperRight,
               std::vector<size_t>&& numElements, std::vector<RangeType>&& values, std::string nm = static_id())
    : lowerLeft_(new std::vector<DomainFieldType>(std::move(lowerLeft)))
    , upperRight_(new std::vector<DomainFieldType>(std::move(upperRight)))
    , numElements_(new std::vector<size_t>(std::move(numElements)))
    , values_(new std::vector<RangeType>(std::move(values)))
    , name_(nm)
  {
    // checks
    if (lowerLeft_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "lowerLeft too small (is " << lowerLeft_->size() << ", should be " << dimDomain << ")");
    if (upperRight_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "upperRight too small (is " << upperRight_->size() << ", should be " << dimDomain << ")");
    if (numElements_->size() < dimDomain)
      DUNE_THROW(Dune::RangeError,
                 "numElements too small (is " << numElements_->size() << ", should be " << dimDomain << ")");
    size_t totalSubdomains = 1;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      const auto& ll = (*lowerLeft_)[dd];
      const auto& ur = (*upperRight_)[dd];
      const auto& ne = (*numElements_)[dd];
      if (!(ll < ur))
        DUNE_THROW(Dune::RangeError, "lowerLeft has to be elementwise smaller than upperRight!");
      totalSubdomains *= ne;
    }
    if (values_->size() < totalSubdomains)
      DUNE_THROW(Dune::RangeError,
                 "values too small (is " << values_->size() << ", should be " << totalSubdomains << ")");
  } // Checkerboard(...)

  Checkerboard(const ThisType& other)
    : lowerLeft_(other.lowerLeft_)
    , upperRight_(other.upperRight_)
    , numElements_(other.numElements_)
    , values_(other.values_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      lowerLeft_   = other.lowerLeft_;
      upperRight_  = other.upperRight_;
      numElements_ = other.numElements_;
      values_      = other.values_;
      name_        = other.name_;
    }
    return *this;
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::static_id() + ".checkerboard";
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    // decide on the subdomain the center of the entity belongs to
    const auto center = entity.geometry().center();
    std::vector<size_t> whichPartition(dimDomain, 0);
    const auto& ll = *lowerLeft_;
    const auto& ur = *upperRight_;
    const auto& ne = *numElements_;
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      // for points that are on upperRight_[d], this selects one partition too much
      // so we need to cap this
      whichPartition[dd] =
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), ne[dd] - 1);
    }
    size_t subdomain = 0;
    if (dimDomain == 1)
      subdomain = whichPartition[0];
    else if (dimDomain == 2)
      subdomain = whichPartition[0] + whichPartition[1] * ne[0];
    else
      subdomain = whichPartition[0] + whichPartition[1] * ne[0] + whichPartition[2] * ne[1] * ne[0];
    // return the component that belongs to the subdomain
    return std::unique_ptr<Localfunction>(new Localfunction(entity, (*values_)[subdomain]));
  } // ... local_function(...)

private:
  std::shared_ptr<const std::vector<DomainFieldType>> lowerLeft_;
  std::shared_ptr<const std::vector<DomainFieldType>> upperRight_;
  std::shared_ptr<const std::vector<size_t>> numElements_;
  std::shared_ptr<const std::vector<RangeType>> values_;
  std::string name_;
}; // class Checkerboard


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(etype, ddim)                                                   \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 1)                                                  \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 2)                                                  \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS(etype, ddim, rdim)                                         \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 1)                                        \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                        \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                              \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                       \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                           \
  DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                     \
  extern template class Dune::Stuff::Functions::Checkerboard<etype, dftype, ddim, rftype, rdim, rcdim>;

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_CHECKERBOARD_LIST_DIMRANGE
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
