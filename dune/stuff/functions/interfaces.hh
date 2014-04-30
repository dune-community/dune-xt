// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

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

#include <dune/stuff/common/memory.hh>

#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/vtk.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

namespace Dune {
namespace Stuff {
namespace Functions {

// forwards, include is below
#if HAVE_DUNE_GRID
template <class GridViewType, int dimRange>
class VisualizationAdapter;
#endif // HAVE_DUNE_GRID

template <class MinuendType, class SubtrahendType>
class Difference;
}

namespace Tags {

class LocalizableFunction
{
};
}


/**
 *  \brief  Interface for a set of globalvalued functions, which can be evaluated locally on one Entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class LocalfunctionSetInterface
{
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

  typedef Functions::Difference<ThisType, ThisType> DifferenceType;

  virtual ~LocalizableFunctionInterface()
  {
  }

  static std::string static_id()
  {
    return "dune.stuff.function";
  }

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "dune.stuff.function";
  }
  /* @} */

  DifferenceType operator-(const ThisType& other) const
  {
    return DifferenceType(*this, other);
  }

#if HAVE_DUNE_GRID
  /**
   * \note  We use the SubsamplingVTKWriter (which is better for higher orders) by default. This means that the grid you
   *        see in the visualization is a refinement of the actual grid!
   */
  template <class GridViewType>
  void visualize(const GridViewType& grid_view, const std::string filename, const bool subsampling = true,
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    if (filename.empty())
      DUNE_THROW(RangeError, "Empty filename given!");
    auto adapter = std::make_shared<Stuff::Functions::VisualizationAdapter<GridViewType, dimRange>>(*this);
    if (subsampling) {
      SubsamplingVTKWriter<GridViewType> vtk_writer(grid_view, VTK::nonconforming);
      vtk_writer.addVertexData(adapter);
      vtk_writer.write(filename, vtk_output_type);
    } else {
      VTKWriter<GridViewType> vtk_writer(grid_view, VTK::nonconforming);
      vtk_writer.addVertexData(adapter);
      vtk_writer.write(filename, vtk_output_type);
    }
  } // ... visualize(...)
#endif // HAVE_DUNE_GRID
}; // class LocalizableFunctionInterface


/**
 * base class for global matrix-valued valued functions that provides automatic local functions via
 * LocalizableFunctionInterface
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class GlobalFunctionInterface
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      LocalfunctionType;
  typedef typename LocalfunctionType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = LocalfunctionType::dimDomain;
  typedef typename LocalfunctionType::DomainType DomainType;

  typedef typename LocalfunctionType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = LocalfunctionType::dimRange;
  static const unsigned int dimRangeCols = LocalfunctionType::dimRangeCols;
  typedef typename LocalfunctionType::RangeType RangeType;

  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  virtual ~GlobalFunctionInterface()
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
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> LocalfunctionType;
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;
#if HAVE_DUNE_FEM
  typedef typename Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                       GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                               rangeDim, 1>>::JacobianRangeType JacobianRangeType;
#else
  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;
#endif

  virtual ~GlobalFunctionInterface()
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

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "You have to imlement it if you intend to use it!");
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
}; // class GlobalFunctionInterface< ..., 1 >


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
    return "dune.stuff.function";
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


} // namespace Stuff
} // namespace Dune

#include "default.hh"
#include "combined.hh"

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
#define DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(etype, ddim)                                                      \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGE(Dune::Stuff::LocalfunctionSetInterface, etype, ddim)                   \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGE(Dune::Stuff::LocalfunctionInterface, etype, ddim)                      \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGE(Dune::Stuff::LocalizableFunctionInterface, etype, ddim)

#define DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGE(cname, etype, ddim)                                              \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGECOLS(cname, etype, ddim, 1)                                             \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGECOLS(cname, etype, ddim, 2)                                             \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGECOLS(cname, etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGECOLS(cname, etype, ddim, rdim)                                    \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 1)                                   \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 2)                                   \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, rcdim)                         \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_RANGEFIELDTYPES(cname, etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_RANGEFIELDTYPES(cname, etype, dftype, ddim, rdim, rcdim)                  \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LAST_EXPANSION(cname, etype, dftype, ddim, double, rdim, rcdim)                      \
  DUNE_STUFF_FUNCTIONS_INTERFACES_LAST_EXPANSION(cname, etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_INTERFACES_LAST_EXPANSION(cname, etype, dftype, ddim, rftype, rdim, rcdim)                \
  extern template class cname<etype, dftype, ddim, rftype, rdim, rcdim>;

#include <dune/stuff/grid/fakeentity.hh>

typedef Dune::Stuff::Grid::FakeEntity<1> DuneStuffFunctionsInterfacesFake1dEntityType;
typedef Dune::Stuff::Grid::FakeEntity<2> DuneStuffFunctionsInterfacesFake2dEntityType;
typedef Dune::Stuff::Grid::FakeEntity<3> DuneStuffFunctionsInterfacesFake3dEntityType;

DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesFake1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesFake2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesFake3dEntityType, 3)

#if HAVE_DUNE_GRID

#include <dune/grid/sgrid.hh>

typedef Dune::SGrid<1, 1>::Codim<0>::Entity DuneStuffFunctionsInterfacesSGrid1dEntityType;
typedef Dune::SGrid<2, 2>::Codim<0>::Entity DuneStuffFunctionsInterfacesSGrid2dEntityType;
typedef Dune::SGrid<3, 3>::Codim<0>::Entity DuneStuffFunctionsInterfacesSGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<1>::Codim<0>::Entity DuneStuffFunctionsInterfacesYaspGrid1dEntityType;
typedef Dune::YaspGrid<2>::Codim<0>::Entity DuneStuffFunctionsInterfacesYaspGrid2dEntityType;
typedef Dune::YaspGrid<3>::Codim<0>::Entity DuneStuffFunctionsInterfacesYaspGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#ifdef ALUGRID_CONFORM
#define DUNE_STUFF_FUNCTIONS_INTERFACES_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#define ALUGRID_CONFORM 1
#endif
#ifdef ENABLE_ALUGRID
#define DUNE_STUFF_FUNCTIONS_INTERFACES_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#define ENABLE_ALUGRID 1
#endif

#include <dune/grid/alugrid.hh>

typedef Dune::ALUSimplexGrid<2, 2>::Codim<0>::Entity DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType;
typedef Dune::ALUSimplexGrid<3, 3>::Codim<0>::Entity DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType;
typedef Dune::ALUCubeGrid<3, 3>::Codim<0>::Entity DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType;

DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#ifdef DUNE_STUFF_FUNCTIONS_INTERFACES_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#undef ALUGRID_CONFORM
#endif
#ifdef DUNE_STUFF_FUNCTIONS_INTERFACES_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#undef ENABLE_ALUGRID
#endif

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_INTERFACES_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_DIMRANGE
#undef DUNE_STUFF_FUNCTIONS_INTERFACES_LIST_CLASSES

extern template class Dune::Stuff::FunctionInterface<double, 1, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 1, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 1, double, 3>;

extern template class Dune::Stuff::FunctionInterface<double, 2, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 2, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 2, double, 3>;

extern template class Dune::Stuff::FunctionInterface<double, 3, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 3, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 3, double, 3>;
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
