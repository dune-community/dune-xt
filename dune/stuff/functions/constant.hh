#ifndef DUNE_STUFF_FUNCTION_CONSTANT_HH
#define DUNE_STUFF_FUNCTION_CONSTANT_HH

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimCols, int rangeDimRows>
class FunctionConstantBase
    : public FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimCols, rangeDimRows>,
      public TimedependentFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimCols, rangeDimRows>
{
  typedef FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimCols, rangeDimRows> BaseType;

public:
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  static const unsigned int dimRangeRows = BaseType::dimRangeRows;
  typedef typename BaseType::RangeType RangeType;

  FunctionConstantBase(const RangeFieldType& constant)
    : constant_(constant)
  {
  }

  FunctionConstantBase(const RangeType& constant)
    : constant_(constant)
  {
  }

  static const std::string id()
  {
    return BaseType::id() + ".constant";
  }

  virtual int order() const
  {
    return 0;
  }

  virtual std::string name() const
  {
    return "function.constant";
  }

  virtual void evaluate(const DomainType& /*arg*/, RangeType& ret) const
  {
    ret = constant_;
  }

  virtual void evaluate(const DomainType& /*arg*/, const double& /*t*/, RangeType& ret) const
  {
    ret = constant_;
  }

private:
  const RangeType constant_;
}; // class FunctionConstantBase


// forward, to allow for specialization
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class FunctionConstant
{
public:
  FunctionConstant() = delete;
}; // class FunctionConstant


template <class DomainFieldImp, int domainDim, class RangeFieldImp>
class FunctionConstant<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
    : public FunctionConstantBase<DomainFieldImp, domainDim, RangeFieldImp, 1, 1>
{
  typedef FunctionConstantBase<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> BaseType;

public:
  typedef FunctionConstant<DomainFieldImp, domainDim, RangeFieldImp, 1, 1> ThisType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::RangeType RangeType;

  FunctionConstant(const RangeFieldType& constant)
    : BaseType(constant)
  {
  }

  FunctionConstant(const RangeType& constant)
    : BaseType(constant)
  {
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

  static ThisType* create(const DSC::ExtendedParameterTree settings)
  {
    return new ThisType(settings.get<RangeFieldType>("value", RangeFieldType(0)));
  } // ... create(...)
}; // class FunctionConstant< ..., 1 >


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_CONSTANT_HH
