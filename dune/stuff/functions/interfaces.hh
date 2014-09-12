// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>
#include <ostream>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/vtk.hh>
#include <dune/stuff/common/filesystem.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif // HAVE_DUNE_FEM

#if HAVE_DUNE_PDELAB
#include <dune/typetree/nodetags.hh>
#include <dune/pdelab/common/function.hh>
#endif // HAVE_DUNE_PDELAB
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/exceptions.hh>


namespace Dune {
namespace Stuff {
namespace Functions {

// forwards, include is below
#if HAVE_DUNE_GRID


template <class GridViewType, int dimRange, int dimRangeCols = 1>
class VisualizationAdapter;


#endif // HAVE_DUNE_GRID


template <class MinuendType, class SubtrahendType>
class Difference;


template <class LeftSummandType, class RightSummandType>
class Sum;


} // namespace Functions
namespace Tags {


class LocalizableFunction
{
};


} // namespace Tags


/**
 *  \brief  Interface for a set of globalvalued functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class LocalfunctionSetInterface
{
  static_assert(EntityImp::dimension == domainDim, "Dimensions do not match!");

  template <class RangeFieldType, int dimRange, int dimRangeCols>
  struct RangeTypeSelector
  {
    typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimRangeCols> type;
  };

  template <class RangeFieldType, int dimRange>
  struct RangeTypeSelector<RangeFieldType, dimRange, 1>
  {
    typedef Dune::FieldVector<RangeFieldType, dimRange> type;
  };

  template <int dimDomain, class RangeFieldType, int dimRange, int dimRangeCols>
  struct JacobianRangeTypeSelector
  {
    typedef double type;
  };

  template <int dimDomain, class RangeFieldType, int dimRange>
  struct JacobianRangeTypeSelector<dimDomain, RangeFieldType, dimRange, 1>
  {
    typedef Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> type;
  };

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;
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
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
                            xx
#else
/*xx*/
#endif
                        ) const
  {
#ifndef DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity().type());
    return reference_element.checkInside(xx);
#else // DUNE_STUFF_FUNCTIONS_DISABLE_CHECKS
    return true;
#endif
  }

  const EntityType& entity_;
}; // class LocalfunctionSetInterface


/**
 *  \brief  Interface for functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class LocalfunctionInterface
    : public LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalfunctionSetInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;

public:
  typedef EntityImp EntityType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
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
  virtual size_t size() const DS_FINAL
  {
    return 1;
  }

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const DS_FINAL
  {
    assert(ret.size() >= 1);
    evaluate(xx, ret[0]);
  }

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const DS_FINAL
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
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class LocalizableFunctionInterface : public IsLocalizableFunction, public Tags::LocalizableFunction
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      ThisType;

public:
  typedef EntityImp EntityType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = rangeDimCols;

  typedef LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols>
      LocalfunctionType;

  typedef typename LocalfunctionType::DomainType DomainType;
  typedef typename LocalfunctionType::RangeType RangeType;
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  static const bool available = false;

  typedef Functions::Difference<ThisType, ThisType> DifferenceType;
  typedef Functions::Sum<ThisType, ThisType> SumType;

  virtual ~LocalizableFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* DUNE_DEPRECATED_MSG("Will be removed, does not work (05.09.2014)!") copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string type() const
  {
    return "stuff.function";
  }

  virtual std::string name() const
  {
    return "stuff.function";
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

#if HAVE_DUNE_GRID
  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default. This means that the grid you
   *        see in the visualization is a refinement of the actual grid!
   */
  template <class GridViewType>
  void visualize(const GridViewType& grid_view, const std::string path, const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (path.empty())
      DUNE_THROW(RangeError, "Empty path given!");
    const auto directory = DSC::directoryOnly(path);
    const auto filename  = DSC::filenameOnly(path);
    auto adapter = std::make_shared<Functions::VisualizationAdapter<GridViewType, dimRange, dimRangeCols>>(*this);
    std::unique_ptr<VTKWriter<GridViewType>> vtk_writer =
        subsampling ? DSC::make_unique<SubsamplingVTKWriter<GridViewType>>(grid_view, VTK::nonconforming)
                    : DSC::make_unique<VTKWriter<GridViewType>>(grid_view, VTK::nonconforming);
    vtk_writer->addVertexData(adapter);
    DSC::testCreateDirectory(directory);
    if (MPIHelper::getCollectiveCommunication().size() == 1)
      vtk_writer->write(path, vtk_output_type);
    else
      vtk_writer->pwrite(filename, directory, "", vtk_output_type);
  } // ... visualize(...)
#endif // HAVE_DUNE_GRID

  virtual void report(std::ostream& out, const std::string prefix = "") const
  {
    out << prefix << "function '" << name() << "' (of type " << type() << ")";
  }

private:
  template <class T>
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*function*/);
}; // class LocalizableFunctionInterface


template <class E, class D, int d, class R, int r, int rC>
std::ostream& operator<<(std::ostream& out, const LocalizableFunctionInterface<E, D, d, R, r, rC>& function)
{
  function.report(out);
  return out;
} // ... operator<<(...)


template <class OtherEntityImp, class GlobalFunctionImp>
struct TransferredGlobalFunction;


/**
 * base class for global matrix-valued valued functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class GlobalFunctionInterface
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunctionInterface()
  {
  }

  virtual ThisType* DUNE_DEPRECATED_MSG("Will be removed, does not work (05.09.2014)!") copy() const
  {
    DUNE_THROW(NotImplemented, "not needed, no meaningful default implementation possible -> exception");
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

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const DS_OVERRIDE DS_FINAL
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

private:
  class Localfunction : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity, const ThisType& global_function)
      : LocalfunctionType(entity)
      , geometry_(entity.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const DS_OVERRIDE DS_FINAL
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
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                 GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim,
                                                         1>>
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
      ,
      public TypeTree::LeafNode,
      public PDELab::
          FunctionInterface<PDELab::FunctionTraits<DomainFieldImp, domainDim, FieldVector<DomainFieldImp, domainDim>,
                                                   RangeFieldImp, rangeDim, FieldVector<RangeFieldImp, rangeDim>>,
                            GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                               rangeDim, 1>>::JacobianRangeType JacobianRangeType;
#else
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
#endif

  virtual ~GlobalFunctionInterface()
  {
  }

  virtual ThisType* DUNE_DEPRECATED_MSG("Will be removed, does not work (05.09.2014)!") copy() const
  {
    DUNE_THROW(NotImplemented, "not needed, no meaningful default implementation possible -> exception");
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

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityImp& entity) const DS_OVERRIDE DS_FINAL
  {
    return Common::make_unique<Localfunction>(entity, *this);
  }

private:
  class Localfunction : public LocalfunctionType
  {
  public:
    Localfunction(const EntityImp& entity, const ThisType& global_function)
      : LocalfunctionType(entity)
      , geometry_(entity.geometry())
      , global_function_(global_function)
    {
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.evaluate(xx_global, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      const auto xx_global = geometry_.global(xx);
      global_function_.jacobian(xx_global, ret);
    }

    virtual size_t order() const DS_OVERRIDE DS_FINAL
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
struct TransferredGlobalFunction
    : public GlobalFunctionInterface<OtherEntityImp, typename GlobalFunctionImp::DomainFieldType,
                                     GlobalFunctionImp::dimDomain, typename GlobalFunctionImp::RangeFieldType,
                                     GlobalFunctionImp::dimRange, GlobalFunctionImp::dimRangeCols>
{

  typedef GlobalFunctionInterface<OtherEntityImp, typename GlobalFunctionImp::DomainFieldType,
                                  GlobalFunctionImp::dimDomain, typename GlobalFunctionImp::RangeFieldType,
                                  GlobalFunctionImp::dimRange, GlobalFunctionImp::dimRangeCols> BaseType;
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

  const GlobalFunctionImp& function_;
};

//! Utility to generate a complete Function Type from an existing one and a template
template <class FunctionImp, template <class, class, int, class, int, int> class OutTemplate>
struct FunctionTypeGenerator
{
  typedef OutTemplate<typename FunctionImp::EntityType, typename FunctionImp::DomainFieldType, FunctionImp::dimDomain,
                      typename FunctionImp::RangeFieldType, FunctionImp::dimRange, FunctionImp::dimRangeCols> type;
};

} // namespace Stuff
} // namespace Dune

#include "default.hh"
#include "combined.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
