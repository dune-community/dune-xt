#ifndef DUNE_STUFF_FUNCTION_FEMADAPTER_HH
#define DUNE_STUFF_FUNCTION_FEMADAPTER_HH


#if HAVE_DUNE_FEM

#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/quadrature/quadrature.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

namespace Dune {
namespace Stuff {

/**
 * throw in a LocalizableFunctionInterface derived class and use this adapter in Dune::Fem classes
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FemFunctionAdapter
    : public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                 FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>>,
      public Dune::Fem::HasLocalFunction
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> InterfaceType;
  typedef FemFunctionAdapter<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

  typedef EntityImp EntityType;
  typedef std::unique_ptr<typename InterfaceType::LocalfunctionType> LocalFunctionPtrType;

  //! hide the virtual lf ptr
  class LocalFunction
  {
  public:
    LocalFunction(LocalFunctionPtrType ptr)
      : lf_ptr_(std::move(ptr))
    {
    }

    void evaluate(const typename InterfaceType::DomainType& xx, typename InterfaceType::RangeType& ret) const
    {
      lf_ptr_->evaluate(xx, ret);
    }

    template <class Quadrature>
    void evaluate(const Dune::Fem::QuadraturePointWrapper<Quadrature>& pw, typename InterfaceType::RangeType& ret) const
    {
      evaluate(pw.quadrature().point(pw.point()), ret);
    }


    void jacobian(const typename InterfaceType::DomainType& xx, typename InterfaceType::JacobianRangeType& ret) const
    {
      lf_ptr_->jacobian(xx, ret);
    }

  private:
    LocalFunctionPtrType lf_ptr_;
  };

public:
  typedef LocalFunction LocalFunctionType;

  FemFunctionAdapter(const InterfaceType& function)
    : function_(function)
  {
  }

  LocalFunction localFunction(const EntityType& entity) const
  {
    return LocalFunction(function_.local_function(entity));
  }

private:
  const InterfaceType& function_;
};

template <class F>
FemFunctionAdapter<typename F::EntityType, typename F::DomainFieldType, F::dimDomain, typename F::RangeFieldType,
                   F::dimRange>
femFunctionAdapter(const F& function)
{
  return FemFunctionAdapter<typename F::EntityType,
                            typename F::DomainFieldType,
                            F::dimDomain,
                            typename F::RangeFieldType,
                            F::dimRange>(function);
}

} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_FEM

#endif // DUNE_STUFF_FUNCTION_FEMADAPTER_HH
