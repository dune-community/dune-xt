#ifndef DUNE_STUFF_FUNCTION_GLOBAL_HH
#define DUNE_STUFF_FUNCTION_GLOBAL_HH

#include <memory>
#include <string>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace Stuff {

/**
 * base class for global valued functions that provides automatic local functions via LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class GlobalFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
#if HAVE_DUNE_FEM
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                 GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>>
#endif // HAVE_DUNE_FEM
{
  typedef GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;
  typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      LocalfunctionBaseType;

public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                      rangeDim>>::JacobianRangeType JacobianRangeType;
#else
  typedef typename LocalfunctionBaseType::JacobianRangeType JacobianRangeType;
#endif

  virtual ~GlobalFunction()
  {
  }

  virtual std::string name() const
  {
    return "";
  }

  virtual ThisType* copy() const
  {
    return new ThisType(*this);
  }

  virtual size_t order() const
  {
    return std::numeric_limits<size_t>::max();
  }

  class Localfunction : public LocalfunctionBaseType
  {
  public:
    Localfunction(const EntityImp& entity, const ThisType& global_function)
      : LocalfunctionBaseType(entity)
      , geometry_(entity.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const
    {
      return global_function_.order();
    }

  private:
    const typename EntityImp::Geometry& geometry_;
    const ThisType& global_function_;
  };

  virtual std::unique_ptr<LocalfunctionBaseType> local_function(const EntityImp& entity) const
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }
}; // class GlobalFunction


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class GlobalConstantFunction : public GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

public:
  explicit GlobalConstantFunction(const RangeType& constant)
    : constant_(constant)
  {
  }

  explicit GlobalConstantFunction(const RangeFieldImp& value)
    : constant_(value)
  {
  }

  virtual size_t order() const final override
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*x*/, RangeType& ret) const final
  {
    ret = constant_;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const final
  {
    ret = JacobianRangeType(0);
  }

private:
  const RangeType constant_;
};

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_GLOBAL_HH
