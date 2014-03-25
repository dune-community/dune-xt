#ifndef DUNE_STUFF_FUNCTION_GLOBAL_HH
#define DUNE_STUFF_FUNCTION_GLOBAL_HH

#include <memory>
#include <string>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#if HAVE_DUNE_PDELAB
#include <dune/typetree/nodetags.hh>
#include <dune/pdelab/common/function.hh>
#endif

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

namespace Dune {
namespace Stuff {


/**
 * base class for global matrix-valued valued functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class GlobalFunction
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;
  typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      LocalfunctionBaseType;

public:
  typedef typename LocalfunctionBaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = LocalfunctionBaseType::dimDomain;
  typedef typename LocalfunctionBaseType::DomainType DomainType;

  typedef typename LocalfunctionBaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = LocalfunctionBaseType::dimRange;
  static const unsigned int dimRangeCols = LocalfunctionBaseType::dimRangeCols;
  typedef typename LocalfunctionBaseType::RangeType RangeType;

  typedef typename LocalfunctionBaseType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunction()
  {
  }

  virtual std::string name() const
  {
    return "dune.stuff.function.global";
  }

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "not needed, no meaningful default implementation possible -> exception");
  }

  virtual size_t order() const
  {
    return std::numeric_limits<size_t>::max();
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "This does not make sense yet for matrix-valued functions!");
  }

  class Localfunction : public LocalfunctionBaseType
  {
  public:
    Localfunction(const EntityImp& entity, const ThisType& global_function)
      : LocalfunctionBaseType(entity)
      //      , geometry_(entity.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const
    {
      const auto xx_global = /*geometry_*/ this->entity().geometry().global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const
    {
      const auto xx_global = /*geometry_*/ this->entity().geometry().global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const
    {
      return global_function_.order();
    }

  private:
    //      const typename EntityImp::Geometry& geometry_;
    const ThisType& global_function_;
  };

  virtual std::unique_ptr<LocalfunctionBaseType> local_function(const EntityImp& entity) const
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }
}; // class GlobalFunction


/**
 * base class for global valued functions that provides automatic local functions via LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
#if HAVE_DUNE_FEM
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                 GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
      ,
      public TypeTree::LeafNode,
      public PDELab::FunctionInterface<PDELab::FunctionTraits<DomainFieldImp, domainDim,
                                                              FieldVector<DomainFieldImp, domainDim>, RangeFieldImp,
                                                              rangeDim, FieldVector<RangeFieldImp, rangeDim>>,
                                       GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif
{
  typedef GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
      LocalfunctionBaseType;

public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                                      1>>::JacobianRangeType JacobianRangeType;
#else
  typedef typename LocalfunctionBaseType::JacobianRangeType JacobianRangeType;
#endif

  virtual ~GlobalFunction()
  {
  }

  virtual std::string name() const
  {
    return "dune.stuff.function.global";
  }

  virtual ThisType* copy() const
  {
    DUNE_THROW(NotImplemented, "not needed, no meaningful default implementation possible -> exception");
  }

  virtual size_t order() const
  {
    return std::numeric_limits<size_t>::max();
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const = 0;
  virtual void jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  class Localfunction : public LocalfunctionBaseType
  {
  public:
    Localfunction(const EntityImp& entity, const ThisType& global_function)
      : LocalfunctionBaseType(entity)
      //      , geometry_(entity.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const
    {
      const auto xx_global = /*geometry_*/ this->entity().geometry().global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const
    {
      const auto xx_global = /*geometry_*/ this->entity().geometry().global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const
    {
      return global_function_.order();
    }

  private:
    //      const typename EntityImp::Geometry& geometry_;
    const ThisType& global_function_;
  };

  virtual std::unique_ptr<LocalfunctionBaseType> local_function(const EntityImp& entity) const
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }
}; // class GlobalFunction< ..., 1 >


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class GlobalConstantFunction
    : public GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
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

  virtual size_t order() const DS_FINAL DS_OVERRIDE
  {
    return 0;
  }

  virtual void evaluate(const DomainType& /*x*/, RangeType& ret) const DS_FINAL
  {
    ret = constant_;
  }

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& ret) const DS_FINAL
  {
    ret *= 0.0;
  }

private:
  const RangeType constant_;
};

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_GLOBAL_HH
