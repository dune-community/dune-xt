#ifndef DUNE_STUFF_FUNCTION_INTERFACE_hh
#define DUNE_STUFF_FUNCTION_INTERFACE_hh

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#if HAVE_EIGEN
#include <Eigen/Core>
#endif // HAVE_EIGEN

#ifdef HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

namespace Dune {
namespace Stuff {
namespace Function {

#ifdef HAVE_DUNE_FEM
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Interface : public Dune::Fem::Function<Dune::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                             Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>>
{
  typedef Dune::Fem::Function<Dune::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                              Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>> BaseType;

public:
  typedef typename BaseType::DomainFieldType DomainFieldType;

  static const int dimDomain = domainDim;

  typedef typename BaseType::DomainType DomainType;

  static const int dimRange = rangeDim;

  typedef typename BaseType::RangeFieldType RangeFieldType;

  typedef typename BaseType::RangeType RangeType;

#if HAVE_EIGEN
  typedef typename Eigen::Matrix<RangeFieldType, Eigen::Dynamic, 1> EigenVectorType;
#endif // HAVE_EIGEN

  virtual void evaluate(const DomainType&, RangeType&) const = 0;

#if HAVE_EIGEN
  virtual void evaluate(const EigenVectorType& arg, EigenVectorType& ret) const
  {
    // process input
    assert(arg.size() == dimDomain);
    DomainType arg_fvector;
    for (int i = 0; i < dimDomain; ++i)
      arg_fvector[i] = arg(i);
    // evaluate
    RangeType ret_fvector;
    evaluate(arg_fvector, ret_fvector);
    // process output
    assert(ret.size() == dimRange);
    for (int i = 0; i < dimRange; ++i)
      ret(i) = ret_fvector[i];
  }

  virtual EigenVectorType evaluate(const EigenVectorType& arg) const
  {
    EigenVectorType ret;
    evaluate(arg, ret);
    return ret;
  }

  virtual void evaluate(const DomainType& arg, EigenVectorType& ret) const
  {
    // evaluate
    RangeType ret_fvector;
    evaluate(arg, ret_fvector);
    // process output
    assert(ret.size() == dimRange);
    for (int i = 0; i < dimRange; ++i)
      ret(i) = ret_fvector[i];
  }

  virtual void evaluate(const EigenVectorType& arg, RangeType& ret) const
  {
    // process input
    assert(arg.size() == dimDomain);
    DomainType arg_fvector;
    for (int i = 0; i < dimDomain; ++i)
      arg_fvector[i] = arg(i);
    // evaluate
    RangeType ret_fvector;
    evaluate(arg_fvector, ret);
  }
#endif // HAVE_EIGEN

  virtual RangeType evaluate(const DomainType& arg) const
  {
    RangeType ret;
    evaluate(arg, ret);
    return ret;
  }
};

#else // HAVE_DUNE_FEM

/**
 *  \todo This should have been replaced by Dune::Function? Investigate further!
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Interface
{
public:
  typedef DomainFieldImp DomainFieldType;

  static const int dimDomain = domainDim;

  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;

  static const int dimRange = rangeDim;

  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

#if HAVE_EIGEN
  typedef typename Eigen::Matrix<RangeFieldType, Eigen::Dynamic, 1> EigenVectorType;
#endif // HAVE_EIGEN

  virtual void evaluate(const DomainType&, RangeType&) const = 0;

  virtual RangeType evaluate(const DomainType& arg) const
  {
    RangeType ret;
    evaluate(arg, ret);
    return ret;
  }

#if HAVE_EIGEN
  virtual void evaluate(const EigenVectorType& arg, EigenVectorType& ret) const
  {
    // process input
    assert(arg.size() == dimDomain);
    DomainType arg_fvector;
    for (int i = 0; i < dimDomain; ++i)
      arg_fvector[i] = arg(i);
    // evaluate
    RangeType ret_fvector;
    evaluate(arg_fvector, ret_fvector);
    // process output
    assert(ret.size() == dimRange);
    for (int i = 0; i < dimRange; ++i)
      ret(i) = ret_fvector[i];
  }

  virtual EigenVectorType evaluate(const EigenVectorType& arg) const
  {
    EigenVectorType ret;
    evaluate(arg, ret);
    return ret;
  }

  virtual void evaluate(const DomainType& arg, EigenVectorType& ret) const
  {
    // evaluate
    RangeType ret_fvector;
    evaluate(arg, ret_fvector);
    // process output
    assert(ret.size() == dimRange);
    for (int i = 0; i < dimRange; ++i)
      ret(i) = ret_fvector[i];
  }

  virtual void evaluate(const EigenVectorType& arg, RangeType& ret) const
  {
    // process input
    assert(arg.size() == dimDomain);
    DomainType arg_fvector;
    for (int i = 0; i < dimDomain; ++i)
      arg_fvector[i] = arg(i);
    // evaluate
    RangeType ret_fvector;
    evaluate(arg_fvector, ret);
  }
#endif // HAVE_EIGEN
}; // class Interface

#endif // HAVE_DUNE_FEM

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_INTERFACE_hh
