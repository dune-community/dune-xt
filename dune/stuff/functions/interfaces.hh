#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>

#include <dune/common/bartonnackmanifcheck.hh>
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

namespace Dune {
namespace Stuff {


// forward, includes are below
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class LocalizedFunction;


/**
 *  \brief  Interface for matrix valued globalvalued functions, which can be evaluated locally on one Entity.
 *
 *          This is the interface for matrixvalued functions, see the specialization for rangeDimCols = 1 for scalar
 *          and vector valued functions.
 */
template <class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class LocalFunctionInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  const EntityType& entity() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().entity());
    return asImp().entity();
  }

  virtual int order() const
  {
    return -1;
  }

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().evaluate(x, ret));
    asImp().evaluate(x, ret);
  }

  RangeType evaluate(const DomainType& x) const
  {
    RangeType ret(0);
    evaluate(x, ret);
    return ret;
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class LocalFunctionInterface


/**
 *  \brief  Interface for scalar and vector valued globalvalued functions, which can be evaluated locally on one Entity.
 */
template <class Traits, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class LocalFunctionInterface<Traits, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::EntityType EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  virtual ~LocalFunctionInterface()
  {
  }

  const EntityType& entity() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().entity());
    return asImp().entity();
  }

  virtual int order() const
  {
    return -1;
  }

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().evaluate(x, ret));
    asImp().evaluate(x, ret);
  }

  RangeType evaluate(const DomainType& x) const
  {
    RangeType ret(0);
    evaluate(x, ret);
    return ret;
  }

  void jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().jacobian(x, ret));
    asImp().jacobian(x, ret);
  }

  JacobianRangeType jacobian(const DomainType& x) const
  {
    JacobianRangeType ret(0);
    jacobian(x, ret);
    return ret;
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class LocalFunctionInterface< ..., 1 >


/**
 *  \brief  Flag for all functions which provide a localFunction(entity) method.
 *
 *  The derived class has to provide a method with the following signature:
\code
template< class EntityType >
ReturnType localFunction(const EntityType& entity) const
{
  ...
}
\endcode
 *  and a struct to provide the ReturnType of this method:
\code
template< class EntityType >
struct LocalFunction
{
  typedef ... Type;
};
\endcode
 *  You can thus obtain the return type and the local function, if derived is an instance of a derived class
 *  of this class with type derived_type:
\code
typedef typename derived_type::template LocalFunction< EntityType >::Type LocalFunctionType;
const LocalFunctionType localfunction = derived.localFunction(entity);
\endcode
 */
class LocalizableFunction
{
};


/**
 * \brief Interface for matrix valued stationary function.
 * \note  See specialization (rangeDimRows = 1) for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
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
    return "dune.stuff.function";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
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
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
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
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                                         dimRangeCols>>::JacobianRangeType JacobianRangeType;
#else
  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;
#endif

  virtual ~FunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "dune.stuff.function";
  }


  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
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

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "You really have to implement this!");
  }

  virtual JacobianRangeType jacobian(const DomainType& x) const
  {
    JacobianRangeType ret;
    jacobian(x, ret);
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

  static std::string static_id()
  {
    return "dune.stuff.timedependentfunction";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
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

  static std::string static_id()
  {
    return "dune.stuff.timedependentfunction";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
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

#include "local.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
