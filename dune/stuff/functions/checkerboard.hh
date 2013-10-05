#ifndef DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
#define DUNE_STUFF_FUNCTION_CHECKERBOARD_HH

#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Checkerboard
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef Checkerboard<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

  class Localfunction : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

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
      : entity_(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual const EntityType& entity() const override
    {
      return entity_;
    }

    virtual size_t order() const override
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const EntityType& entity_;
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

  virtual ThisType* copy() const override
  {
    return new ThisType(*this);
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
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
          std::min(size_t(std::floor(ne[dd] * ((center[dd] - ll[dd]) / (ur[dd] - ll[dd])))), (ne[dd] - 1));
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


//// forward, to allow for specialization
// template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
// class FunctionCheckerboard
//{
// public:
//  FunctionCheckerboard() = delete;
//};


// template< class DomainFieldImp, int domainDim, class RangeFieldImp >
// class FunctionCheckerboard< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//  : public FunctionCheckerboardBase< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//{
//  typedef FunctionCheckerboardBase< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > BaseType;
// public:
//  typedef FunctionCheckerboard< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 > ThisType;

//  typedef typename BaseType::DomainFieldType  DomainFieldType;
//  static const unsigned int                   dimDomain =BaseType::dimDomain;
//  typedef typename BaseType::DomainType       DomainType;
//  typedef typename BaseType::RangeFieldType   RangeFieldType;
//  typedef typename BaseType::RangeType        RangeType;

//  using BaseType::static_id;

//  FunctionCheckerboard(const DomainType _lowerLeft,
//                       const DomainType _upperRight,
//                       const std::vector< size_t > _numElements,
//                       const std::vector< RangeType > _values,
//                       const std::string _name = static_id())
//    : BaseType(_lowerLeft, _upperRight, _numElements, _values, _name)
//  {}

//  static Dune::ParameterTree defaultSettings(const std::string subName = "")
//  {
//    Dune::ParameterTree description;
//    description["lowerLeft"] = "[0.0; 0.0; 0.0]";
//    description["upperRight"] = "[1.0; 1.0; 1.0]";
//    description["numElements"] = "[2; 2; 2]";
//    description["values"] = "[1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0]";
//    description["name"] = static_id();
//    if (subName.empty())
//      return description;
//    else {
//      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
//      extendedDescription.add(description, subName);
//      return extendedDescription;
//    }
//  } // ... defaultSettings(...)

//  static ThisType* create(const DSC::ExtendedParameterTree settings = defaultSettings())
//  {
//    // get data
//    const std::vector< DomainFieldType > lowerLefts = settings.getVector("lowerLeft", DomainFieldType(0), dimDomain);
//    const std::vector< DomainFieldType > upperRights = settings.getVector("upperRight",
//                                                                             DomainFieldType(1),
//                                                                             dimDomain);
//    const std::vector< size_t > numElements = settings.getVector("numElements", size_t(1), dimDomain);
//    size_t subdomains = 1;
//    for (size_t ii = 0; ii < numElements.size(); ++ii)
//      subdomains *= numElements[ii];
//    const std::vector< RangeFieldType > rangeFieldTypeValues = settings.getVector("values", RangeFieldType(1),
//    subdomains);
//    std::vector< RangeType > rangeTypeValues(rangeFieldTypeValues.size(), RangeType(0));
//    for (size_t ii = 0; ii < rangeFieldTypeValues.size(); ++ii)
//      rangeTypeValues[ii] = rangeFieldTypeValues[ii];
//    // convert and leave the checks to the base constructor
//    DomainType lowerLeft;
//    DomainType upperRight;
//    for (size_t dd = 0; dd < dimDomain; ++dd) {
//      lowerLeft[dd] = lowerLefts[dd];
//      upperRight[dd] = upperRights[dd];
//    }
//    // create and return
//    return new ThisType(lowerLeft, upperRight, numElements, rangeTypeValues);
//  } // ... create(...)
//};


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CHECKERBOARD_HH
