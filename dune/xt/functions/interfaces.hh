// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2015)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014 - 2015)

#ifndef DUNE_XT_FUNCTIONS_INTERFACE_HH
#define DUNE_XT_FUNCTIONS_INTERFACE_HH

#include <memory>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/version.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/vtk.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/common/function.hh>

#include <dune/typetree/nodetags.hh>
#endif

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {

template <class F>
struct is_localizable_function_helper
{
  DXTC_has_typedef_initialize_once(EntityType);
  DXTC_has_typedef_initialize_once(DomainFieldType);
  DXTC_has_typedef_initialize_once(RangeFieldType);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);

  static const bool is_candidate =
      DXTC_has_typedef(EntityType)<F>::value && DXTC_has_typedef(DomainFieldType)<F>::value
      && DXTC_has_typedef(RangeFieldType)<F>::value && DXTC_has_static_member(dimDomain)<F>::value
      && DXTC_has_static_member(dimRange)<F>::value && DXTC_has_static_member(dimRangeCols)<F>::value;
}; // class is_localizable_function_helper

} // namespace internal

// forwards, includes are below
template <class F, bool candidate = internal::is_localizable_function_helper<F>::is_candidate>
struct is_localizable_function;

namespace internal {

// additional argument for member functions to differentiate between dimRangeCols = 1 and dimRangeCols > 1 by
// overloading
template <size_t rangeDimCols>
struct ChooseVariant
{
};

} // namespace internal


template <class GridViewType, size_t dimRange, size_t dimRangeCols = 1>
class VisualizationAdapterFunction;


template <class MinuendType, class SubtrahendType>
class DifferenceFunction;

template <class LeftSummandType, class RightSummandType>
class SumFunction;

template <class LeftSummandType, class RightSummandType>
class ProductFunction;

template <class FunctionImp>
class DivergenceFunction;


/**
 *  \brief Interface for a set of globalvalued functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class LocalfunctionSetInterface
{
  static_assert(EntityImp::dimension == domainDim, "Dimensions do not match!");

  template <class RangeFieldType, size_t dimRange, size_t dimRangeCols>
  struct RangeTypeSelector
  {
    typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimRangeCols> type;
  };

  template <class RangeFieldType, size_t dimRange>
  struct RangeTypeSelector<RangeFieldType, dimRange, 1>
  {
    typedef Dune::FieldVector<RangeFieldType, dimRange> type;
  };

  template <size_t dimDomain, class RangeFieldType, size_t dimRange, size_t dimRangeCols>
  struct JacobianRangeTypeSelector
  {
    typedef Dune::FieldVector<Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain>, dimRangeCols> type;
  };

  template <size_t dimDomain, class RangeFieldType, size_t dimRange>
  struct JacobianRangeTypeSelector<dimDomain, RangeFieldType, dimRange, 1>
  {
    typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> type;
  };

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const constexpr size_t dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const constexpr size_t dimRange = rangeDim;
  static const constexpr size_t dimRangeCols = rangeDimCols;
  typedef typename RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type RangeType;
  typedef typename JacobianRangeTypeSelector<dimDomain, RangeFieldType, dimRange, dimRangeCols>::type JacobianRangeType;

  LocalfunctionSetInterface(const EntityType& ent)
    : entity_(ent)
  {
  }

  virtual ~LocalfunctionSetInterface()
  {
  }

  virtual const EntityType& entity() const
  {
    return entity_;
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
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

protected:
  bool is_a_valid_point(const DomainType&
#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
                            xx
#else
/*xx*/
#endif
                        ) const
  {
#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity().type());
    return reference_element.checkInside(xx);
#else // DUNE_XT_FUNCTIONS_DISABLE_CHECKS
    return true;
#endif
  }

  const EntityType& entity_;
}; // class LocalfunctionSetInterface

/**
 *  \brief  Interface for functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class LocalfunctionInterface
    : public LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;

public:
  typedef EntityImp EntityType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const constexpr size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  LocalfunctionInterface(const EntityType& ent)
    : BaseType(ent)
  {
  }

  virtual ~LocalfunctionInterface()
  {
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented in addition to the ones required from the BaseType.''
   * @{
   **/
  virtual void evaluate(const DomainType& /*xx*/, RangeType& /*ret*/) const = 0;

  virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const = 0;
  /* @} */

  /**
   * \defgroup providedbase ´´These methods are provided by the interface to please LocalfunctionSetInterface.''
   * @{
   **/
  virtual size_t size() const override final
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const override final
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const override final
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

  //! evaluate at N quadrature points into vector of size >= N
  void evaluate(const Dune::QuadratureRule<DomainFieldType, dimDomain>& quadrature, std::vector<RangeType>& ret)
  {
    assert(ret.size() >= quadrature.size());
    std::size_t i = 0;
    for (const auto& point : quadrature)
      evaluate(point.position(), ret[i++]);
  }

  //! jacobian at N quadrature points into vector of size >= N
  void jacobian(const Dune::QuadratureRule<DomainFieldType, dimDomain>& quadrature, std::vector<JacobianRangeType>& ret)
  {
    assert(ret.size() >= quadrature.size());
    std::size_t i = 0;
    for (const auto& point : quadrature)
      jacobian(point.position(), ret[i++]);
  }
  /* @} */
}; // class LocalfunctionInterface

class IsLocalizableFunction
{
};

/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class LocalizableFunctionInterface : public IsLocalizableFunction
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ThisType;

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const constexpr size_t dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const constexpr size_t dimRange = rangeDim;
  static const constexpr size_t dimRangeCols = rangeDimCols;

  typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  static const bool available = false;

  typedef Functions::DifferenceFunction<ThisType, ThisType> DifferenceType;
  typedef Functions::SumFunction<ThisType, ThisType> SumType;
  typedef Functions::DivergenceFunction<ThisType> DivergenceType;

  virtual ~LocalizableFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string type() const
  {
    return "function";
  }

  virtual std::string name() const
  {
    return "function";
  }
  /* @} */

  DifferenceType operator-(const ThisType& other) const
  {
    return DifferenceType(*this, other);
  }

  SumType operator+(const ThisType& other) const
  {
    return SumType(*this, other);
  }

  template <class OtherType>
  typename std::enable_if<is_localizable_function<OtherType>::value,
                          Functions::ProductFunction<ThisType, OtherType>>::type
  operator*(const OtherType& other) const
  {
    return Functions::ProductFunction<ThisType, OtherType>(*this, other);
  }

  DivergenceType divergence() const
  {
    return DivergenceType(*this);
  }

  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default. This means that the grid you
   *        see in the visualization is a refinement of the actual grid!
   */
  template <class GridViewType>
  void visualize(const GridViewType& grid_view,
                 const std::string path,
                 const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (path.empty())
      DUNE_THROW(RangeError, "Empty path given!");
    const auto directory = Common::directory_only(path);
    const auto filename = Common::filename_only(path);
    const auto adapter = std::make_shared<VisualizationAdapterFunction<GridViewType, dimRange, dimRangeCols>>(*this);
    std::unique_ptr<VTKWriter<GridViewType>> vtk_writer =
        subsampling ? Common::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, VTK::nonconforming)
                    : Common::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
    vtk_writer->addVertexData(adapter);
    Common::test_create_directory(directory);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer->write(path, vtk_output_type);
    else
      vtk_writer->pwrite(filename, directory, "", vtk_output_type);
  } // ... visualize(...)

  virtual void report(std::ostream& out, const std::string prefix = "") const
  {
    out << prefix << "function '" << name() << "' (of type " << type() << ")";
  }

private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*function*/);
}; // class LocalizableFunctionInterface

template <class E, class D, size_t d, class R, size_t r, size_t rC>
std::ostream& operator<<(std::ostream& out, const LocalizableFunctionInterface<E, D, d, R, r, rC>& function)
{
  function.report(out);
  return out;
} // ... operator<<(...)

template <class OtherEntityImp, class GlobalFunctionImp>
class TransferredGlobalFunction;

/**
 * base class for global matrix-valued valued functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class GlobalFunctionInterface
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunctionInterface()
  {
  }

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& xx, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "This does not make sense yet for matrix-valued functions!");
  }

  virtual RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret;
    evaluate(xx, ret);
    return ret;
  }

  virtual JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret);
    return ret;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

  virtual std::string type() const override
  {
    return "globalfunction";
  }

  virtual std::string name() const override
  {
    return "globalfunction";
  }

private:
  class Localfunction : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity_in, const ThisType& global_function)
      : LocalfunctionType(entity_in)
      , geometry_(entity_in.geometry())
      , global_function_(global_function)
    {
    }

    virtual ~Localfunction()
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const override final
    {
      return global_function_.order();
    }

  private:
    const typename EntityImp::Geometry geometry_;
    const ThisType& global_function_;
  }; // class Localfunction

public:
  template <class OtherEntityImp>
  struct Transfer
  {
    typedef TransferredGlobalFunction<OtherEntityImp, ThisType> Type;
  };

  template <class OtherEntityImp>
  typename Transfer<OtherEntityImp>::Type transfer() const
  {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFunctionInterface

/**
 * base class for global valued functions that provides automatic local functions via LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
#if HAVE_DUNE_FEM
      ,
      public Dune::Fem::
          Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                   GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
      ,
      public TypeTree::LeafNode,
      public PDELab::
          FunctionInterface<PDELab::FunctionTraits<DomainFieldImp,
                                                   domainDim,
                                                   FieldVector<DomainFieldImp, domainDim>,
                                                   RangeFieldImp,
                                                   rangeDim,
                                                   FieldVector<RangeFieldImp, rangeDim>>,
                            GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::EntityType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::
      Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
               GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>::
          JacobianRangeType JacobianRangeType;
#else
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
#endif

  static_assert(std::is_same<typename LocalfunctionType::RangeType, RangeType>::value, "RangeType mismatch");
  static_assert(std::is_same<typename LocalfunctionType::DomainType, DomainType>::value, "DomainType mismatch");
  static_assert(std::is_same<Dune::FieldVector<DomainFieldImp, domainDim>, DomainType>::value, "DomainType mismatch");

  virtual ~GlobalFunctionInterface()
  {
  }

  virtual size_t order() const = 0;

  virtual void evaluate(const DomainType& x, RangeType& ret) const = 0;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "You have to imlement it if you intend to use it!");
  }

  virtual RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret;
    evaluate(xx, ret);
    return ret;
  }

  virtual JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret);
    return ret;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const override final
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

  virtual std::string type() const override
  {
    return "globalfunction";
  }

  virtual std::string name() const override
  {
    return "globalfunction";
  }

private:
  class Localfunction : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity_in, const ThisType& global_function)
      : LocalfunctionType(entity_in)
      , geometry_(entity_in.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override final
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const override final
    {
      return global_function_.order();
    }

  private:
    const typename EntityImp::Geometry geometry_;
    const ThisType& global_function_;
  }; // class Localfunction

public:
  template <class OtherEntityImp>
  struct Transfer
  {
    typedef TransferredGlobalFunction<OtherEntityImp, ThisType> Type;
  };

  template <class OtherEntityImp>
  typename Transfer<OtherEntityImp>::Type transfer() const
  {
    return typename Transfer<OtherEntityImp>::Type(*this);
  }
}; // class GlobalFunctionInterface< ..., 1 >

template <class OtherEntityImp, class GlobalFunctionImp>
class TransferredGlobalFunction : public GlobalFunctionInterface<OtherEntityImp,
                                                                 typename GlobalFunctionImp::DomainFieldType,
                                                                 GlobalFunctionImp::dimDomain,
                                                                 typename GlobalFunctionImp::RangeFieldType,
                                                                 GlobalFunctionImp::dimRange,
                                                                 GlobalFunctionImp::dimRangeCols>
{
  typedef GlobalFunctionInterface<OtherEntityImp,
                                  typename GlobalFunctionImp::DomainFieldType,
                                  GlobalFunctionImp::dimDomain,
                                  typename GlobalFunctionImp::RangeFieldType,
                                  GlobalFunctionImp::dimRange,
                                  GlobalFunctionImp::dimRangeCols>
      BaseType;

public:
  TransferredGlobalFunction(const GlobalFunctionImp& function)
    : function_(function)
  {
  }

  virtual size_t order() const
  {
    return function_.order();
  }

  virtual void evaluate(const typename BaseType::DomainType& x, typename BaseType::RangeType& ret) const
  {
    function_.evaluate(x, ret);
  }

private:
  const GlobalFunctionImp& function_;
}; // class TransferredGlobalFunction


template <class TimeIndependentFunctionImp, class TimeFieldImp = double>
class TimeDependentFunctionInterface
{
public:
  typedef TimeIndependentFunctionImp TimeIndependentFunctionType;
  typedef TimeFieldImp TimeFieldType;

  static const bool available = false;

  virtual ~TimeDependentFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "stuff.timedependentfunction";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr<TimeIndependentFunctionType> evaluate_at_time(const TimeFieldType t) const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string type() const
  {
    return "stuff.timedependentfunction";
  }

  virtual std::string name() const
  {
    return "stuff.timedependentfunction";
  }
  /* @} */
};

//! Utility to generate a complete Function Type from an existing one and a template
template <class FunctionImp, template <class, class, size_t, class, size_t, size_t> class OutTemplate>
struct FunctionTypeGenerator
{
  typedef OutTemplate<typename FunctionImp::EntityType,
                      typename FunctionImp::DomainFieldType,
                      FunctionImp::dimDomain,
                      typename FunctionImp::RangeFieldType,
                      FunctionImp::dimRange,
                      FunctionImp::dimRangeCols>
      type;
};

template <class F>
struct is_localizable_function<F, true>
    : public std::is_base_of<LocalizableFunctionInterface<typename F::EntityType,
                                                          typename F::DomainFieldType,
                                                          F::dimDomain,
                                                          typename F::RangeFieldType,
                                                          F::dimRange,
                                                          F::dimRangeCols>,
                             F>
{
};

template <class F>
struct is_localizable_function<F, false> : public std::false_type
{
};

} // namespace Functions
} // namespace XT
} // namespace Dune

#include "combined.hh"
#include "default.hh"
#include "derived.hh"

#endif // DUNE_XT_FUNCTIONS_INTERFACE_HH
