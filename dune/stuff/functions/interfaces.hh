#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

namespace Dune {
namespace Stuff {
namespace Function {

// forward, include is below
template <class GridViewType, int dimRange>
class VisualizationAdapter;
}


/**
 *  \brief  Interface for a set of matrix valued globalvalued functions, which can be evaluated locally on one Entity.
 *
 *  \note   see specialization for rangeDimCols = 1 for vector and scalar valued localfunction sets.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class LocalfunctionSetInterface
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  //  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType; // <- this is yet unclear

  virtual ~LocalfunctionSetInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual const EntityType& entity() const = 0;

  virtual size_t size() const = 0;

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::vector<RangeType>& /*ret*/) const = 0;

  //  virtual void jacobian(const DomainType& /*xx*/, std::vector< JacobianRangeType >& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  std::vector<RangeType> evaluate(const DomainType& xx) const
  {
    std::vector<RangeType> ret(size(), RangeType(0));
    evaluate(xx, ret);
    return ret;
  }

  //  std::vector< JacobianRangeType > jacobian(const DomainType& xx) const
  //  {
  //    std::vector< JacobianRangeType > ret(size(), JacobianRangeType(0));
  //    jacobian(xx, ret);
  //    return ret;
  //  }
  /* @} */
}; // class LocalfunctionSetInterface


/**
 *  \brief  Interface for a set of scalar or vector globalvalued functions, which can be evaluated locally on one
 * Entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  virtual ~LocalfunctionSetInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual const EntityType& entity() const = 0;

  virtual size_t size() const = 0;

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& /*xx*/, std::vector<RangeType>& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, std::vector<JacobianRangeType>& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  std::vector<RangeType> evaluate(const DomainType& xx) const
  {
    std::vector<RangeType> ret(size(), RangeType(0));
    evaluate(xx, ret);
    return ret;
  }

  std::vector<JacobianRangeType> jacobian(const DomainType& xx) const
  {
    std::vector<JacobianRangeType> ret(size(), JacobianRangeType(0));
    jacobian(xx, ret);
    return ret;
  }
  /* @} */
}; // class LocalfunctionSetInterface< ...., 1 >


/**
 *  \brief  Interface for matrix valued globalvalued functions, which can be evaluated locally on one Entity.
 *
 *          This is the interface for matrixvalued functions, see the specialization for rangeDimCols = 1 for scalar
 *          and vector valued functions.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class LocalfunctionInterface
    : public LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  //  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType; // <- this is yet unclear

  virtual ~LocalfunctionInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to entity() and order().''
   * @{
   **/
  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/) const = 0;

  //  virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup providedbase ´´These methods are provided by the interface to please LocalfunctionSetInterface.''
   * @{
   **/
  virtual size_t size() const
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  //  virtual void jacobian(const DomainType& xx, std::vector< JacobianRangeType >& ret) const
  //  {
  //    assert(ret.size() >= 1);
  //    jacobian(xx, ret[0]);
  //  }
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret(0);
    evaluate(xx, ret);
    return ret;
  }

  //  JacobianRangeType jacobian(const DomainType& xx) const
  //  {
  //    JacobianRangeType ret(0);
  //    jacobian(xx, ret);
  //    return ret;
  //  }
  /* @} */
}; // class LocalfunctionInterface


/**
 *  \brief  Interface for scalar and vector valued globalvalued functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> JacobianRangeType;

  virtual ~LocalfunctionInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to entity() and and order().''
   * @{
   **/
  virtual void evaluate(const DomainType& /*xx*/, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& ret) const = 0;
  /* @} */

  /**
   * \defgroup providedbase ´´These methods are provided by the interface to please LocalfunctionSetInterface.''
   * @{
   **/
  virtual size_t size() const
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const
  {
    assert(ret.size() >= 1);
    jacobian(xx, ret[0]);
  }
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret(0);
    evaluate(xx, ret);
    return ret;
  }

  JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret(0);
    jacobian(xx, ret);
    return ret;
  }
}; // class LocalfunctionInterface< ..., 1 >


class IsLocalizableFunction
{
};


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class LocalizableFunctionInterface : public IsLocalizableFunction
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      ThisType;

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;

  typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
      LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  //  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  virtual ~LocalizableFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "dune.stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to the ones required by the base classes.''
   * @{
   **/
  virtual std::shared_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
  }
  /* @} */
}; // class LocalizableFunctionInterface


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public IsLocalizableFunction
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = rangeDim;
  static const unsigned int dimRangeCols = 1;

  typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1> LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  virtual ~LocalizableFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "dune.stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to the ones required by the base classes.''
   * @{
   **/
  virtual std::shared_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return static_id();
  }
  /* @} */

  template <class GridViewType>
  void visualize(const GridViewType& grid_view, const std::string filename) const
  {
    if (filename.empty())
      DUNE_THROW(IOError, "Empty filename given!");
    auto adapter = std::make_shared<Function::VisualizationAdapter<GridViewType, dimRange>>(*this);
    VTKWriter<GridViewType> vtk_writer(grid_view, VTK::nonconforming);
    vtk_writer.addVertexData(adapter);
    vtk_writer.write(filename);
  } // ... visualize(...)
}; // class LocalizableFunctionInterface


/**
 * \brief Interface for scalar and vector valued stationary function.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class DUNE_DEPRECATED_MSG("Please derive your functions from LocalizableFunctionInterface in the future!")
    FunctionInterface
#if HAVE_DUNE_FEM
    : public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                 FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>>
#endif // HAVE_DUNE_FEM
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp,
                                                         rangeDim>>::JacobianRangeType JacobianRangeType;
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
}; // class FunctionInterface


/**
 * \brief Interface for scalar and vector valued timedependent functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class TimedependentFunctionInterface
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange = rangeDim;
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
  virtual void evaluate(const DomainType& /*xx*/, const double& /*tt*/, RangeType& /*ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


//! use this to throw a stationary function into an algorithm that expects an instationary one
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows>
struct TimeFunctionAdapter
    : public Dune::Stuff::TimedependentFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows>
{
  typedef Dune::Stuff::FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows> WrappedType;

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

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows>
TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows> timefunctionAdapted(
    const Dune::Stuff::FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows>& wrapped)
{
  return TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows>(wrapped);
}


} // namespace Stuff
} // namespace Dune

#include "visualization.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
