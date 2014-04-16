#ifndef DUNE_STUFF_FUNCTIONS_CONSTANT_HH
#define DUNE_STUFF_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Constant
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

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

    Localfunction(const EntityType& ent, const std::shared_ptr<const RangeType> value)
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

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret = *value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const std::shared_ptr<const RangeType> value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

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

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, value_));
  }

private:
  std::shared_ptr<const RangeType> value_;
  std::string name_;
}; // class Constant


} // namespace Function
} // namespace Stuff
} // namespace Dune

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(etype, ddim)                                                       \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 1)                                                      \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 2)                                                      \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS(etype, ddim, rdim)                                             \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 1)                                            \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                            \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                  \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                           \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                               \
  DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                         \
  extern template class Dune::Stuff::Function::Constant<etype, dftype, ddim, rftype, rdim, rcdim>;

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesFake3dEntityType, 3)

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_CONSTANT_LIST_DIMRANGE
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTIONS_CONSTANT_HH
