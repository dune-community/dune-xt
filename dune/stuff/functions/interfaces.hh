#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#include <dune/stuff/common/color.hh>
#include <dune/stuff/fem/namespace.hh>

#include <dune/stuff/localfunction/interface.hh>

#include "local.hh"

namespace Dune {
namespace Stuff {


/**
 * \brief Interface for matrix valued stationary function.
 * \note  See specialization (rangeDimRows = 1) for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class FunctionInterface : public LocalizableFunction
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  virtual ~FunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "function";
  }

  virtual std::string id() const
  {
    return static_id();
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return id();
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*x*/, RangeType& /*ret*/) const = 0;
  /* @} */

  virtual RangeType evaluate(const DomainType& x) const
  {
    RangeType ret;
    evaluate(x, ret);
    return ret;
  }

  template <class EntityType>
  struct LocalFunction
  {
    typedef LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols> Type;
  };

  template <class EntityType>
  LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
  localFunction(const EntityType& entity) const
  {
    return LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>(
        *this, entity);
  }
}; // class FunctionInterface


/**
 * \brief Interface for scalar and vector valued stationary function.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunction
#if HAVE_DUNE_FEM
#if DUNE_FEM_IS_LOCALFUNCTIONS_COMPATIBLE
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
#else
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
#endif
                                 FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif // HAVE_DUNE_FEM
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  virtual ~FunctionInterface()
  {
  }

  static std::string id()
  {
    return "function";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return id();
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*x*/, RangeType& /*ret*/) const = 0;
  /* @} */

  virtual RangeType evaluate(const DomainType& x) const
  {
    RangeType ret;
    evaluate(x, ret);
    return ret;
  }

  template <class EntityType>
  struct LocalFunction
  {
    typedef LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols> Type;
  };

  template <class EntityType>
  LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
  localFunction(const EntityType& entity) const
  {
    return LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>(
        *this, entity);
  }
}; // class FunctionInterface


/**
 * \brief Interface for matrix valued timedependent function.
 * \note  See spezialization (rangeDimCols = 1) for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class TimedependentFunctionInterface
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  virtual ~TimedependentFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.timedependent";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "function.timedependent";
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const double& /*_t*/, RangeType& /*_ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


/**
 * \brief Interface for scalar and vector valued timedependent functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class TimedependentFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  virtual ~TimedependentFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.timedependent";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "function.timedependent";
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const double& /*_t*/, RangeType& /*_ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


//! use this to throw a stationary function into an algorithm that expects an instationary one
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
struct TimeFunctionAdapter
    : public Dune::Stuff::TimedependentFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows,
                                                         rangeDimCols>
{
  typedef Dune::Stuff::FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      WrappedType;

  TimeFunctionAdapter(const WrappedType& wr)
    : wrapped_(wr)
  {
  }

  virtual void evaluate(const typename WrappedType::DomainType& x, typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  virtual void evaluate(const typename WrappedType::DomainType& x, const double& /*t*/,
                        typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  const WrappedType& wrapped_;
};

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> timefunctionAdapted(
    const Dune::Stuff::FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>& wrapped)
{
  return TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>(wrapped);
}


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
