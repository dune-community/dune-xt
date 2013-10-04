#ifndef DUNE_STUFF_FUNCTIONS_CONSTANT_HH
#define DUNE_STUFF_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


/**
 *  \note Only implemented scalar and vector valued at the moment!
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Constant : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

  class Localfunction : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange = BaseType::dimRange;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const std::shared_ptr<const RangeType> value)
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
      ret = *value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const EntityType& entity_;
    const std::shared_ptr<const RangeType> value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = BaseType::dimRange;
  typedef typename LocalfunctionType::RangeType RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".constant";
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["value"] = "1.0";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const DSC::ExtendedParameterTree settings = defaultSettings())
  {
    return new ThisType(settings.get<RangeFieldType>("value", RangeFieldType(0)));
  } // ... create(...)

  Constant(const RangeFieldType& val, const std::string nm = static_id())
    : value_(std::make_shared<RangeType>(val))
    , name_(nm)
  {
  }

  Constant(const RangeType& val, const std::string nm = static_id())
    : value_(std::make_shared<RangeType>(val))
    , name_(nm)
  {
  }

  Constant(const ThisType& other)
    : value_(other.value_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      value_ = other.value_;
      name_  = other.name_;
    }
    return *this;
  }

private:
  Constant(const std::shared_ptr<const RangeType>& val, const std::string nm = static_id())
    : value_(val)
    , name_(nm)
  {
  }

public:
  virtual ThisType* copy() const override
  {
    return new ThisType(value_, name_);
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::shared_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return std::shared_ptr<Localfunction>(new Localfunction(entity, value_));
  }

private:
  std::shared_ptr<const RangeType> value_;
  std::string name_;
}; // class ConstantBase


//// forward, to allow for specialization
// template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1 >
// class FunctionConstant
//  : public FunctionConstantBase< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >
//{
//  typedef FunctionConstantBase< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >  BaseType;
// public:
//  typedef FunctionConstant< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >      ThisType;

//  typedef typename BaseType::RangeFieldType RangeFieldType;
//  typedef typename BaseType::RangeType      RangeType;

//  FunctionConstant(const RangeFieldType& constant)
//    : BaseType(constant)
//  {}

//  FunctionConstant(const RangeType& constant)
//    : BaseType(constant)
//  {}

//  using BaseType::localFunction;
//}; // class FunctionConstant


// template< class DomainFieldImp, int domainDim, class RangeFieldImp >
// class FunctionConstant< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//  : public FunctionConstantBase< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >
//{
//  typedef FunctionConstantBase< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >  BaseType;
// public:
//  typedef FunctionConstant< DomainFieldImp, domainDim, RangeFieldImp, 1, 1 >      ThisType;

//  typedef typename BaseType::RangeFieldType RangeFieldType;
//  typedef typename BaseType::RangeType      RangeType;

//  FunctionConstant(const RangeFieldType& constant)
//    : BaseType(constant)
//  {}

//  FunctionConstant(const RangeType& constant)
//    : BaseType(constant)
//  {}

//  using BaseType::static_id;

//  static Dune::ParameterTree defaultSettings(const std::string subName = "")
//  {
//    Dune::ParameterTree description;
//    description["value"] = "1.0";
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
//    return new ThisType(settings.get< RangeFieldType >("value", RangeFieldType(0)));
//  } // ... create(...)
//}; // class FunctionConstant< ..., 1 >


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_CONSTANT_HH
