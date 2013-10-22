// This file is part of the dune-stuff project:
//   http://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
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

#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

namespace Dune {
namespace Stuff {
#if HAVE_DUNE_GRID
namespace Function {

// forward, include is below
template <class GridViewType, int dimRange>
class VisualizationAdapter;
}
#endif // HAVE_DUNE_GRID


/**
 *  \brief  Interface for a set of matrix valued functions, which can be evaluated locally on one Entity.
 *
 *  \note   see specialization for rangeDimCols = 1 for vector and scalar valued localfunction sets.
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

  LocalfunctionSetInterface(const EntityType& ent);

  virtual ~LocalfunctionSetInterface();

  virtual const EntityType& entity() const;

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
  std::vector<RangeType> evaluate(const DomainType& xx) const;

  std::vector<JacobianRangeType> jacobian(const DomainType& xx) const;
  /* @} */

protected:
  bool is_a_valid_point(const DomainType& xx) const;

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

  LocalfunctionInterface(const EntityType& ent);

  virtual ~LocalfunctionInterface();

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
  virtual size_t size() const;

  virtual void evaluate(const DomainType& xx, std::vector<RangeType>& ret) const;

  virtual void jacobian(const DomainType& xx, std::vector<JacobianRangeType>& ret) const;
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface.''
   * @{
   **/
  RangeType evaluate(const DomainType& xx) const;

  JacobianRangeType jacobian(const DomainType& xx) const;
  /* @} */
}; // class LocalfunctionInterface


class IsLocalizableFunction
{
};


/**
 * \brief Interface for functions which provide a LocalfunctionInterface for an entity.
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class LocalizableFunctionInterface : public IsLocalizableFunction
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

  virtual ~LocalizableFunctionInterface();

  static std::string static_id();

  /**
   * \defgroup haveto ´´These methods have to be implemented.''
   * @{
   **/
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& /*entity*/) const = 0;

  virtual ThisType* copy() const = 0;
  /* @} */

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const;
/* @} */

#if HAVE_DUNE_GRID
  template <class GridViewType>
  void visualize(const GridViewType& grid_view, const std::string filename) const
  {
    if (filename.empty())
      DUNE_THROW(RangeError, "Empty filename given!");
    auto adapter = std::make_shared<Stuff::Function::VisualizationAdapter<GridViewType, dimRange>>(*this);
    VTKWriter<GridViewType> vtk_writer(grid_view, VTK::nonconforming);
    vtk_writer.addVertexData(adapter);
    vtk_writer.write(filename);
  } // ... visualize(...)
#endif // HAVE_DUNE_GRID
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

  virtual ~FunctionInterface();

  static std::string static_id();

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const;

  virtual int order() const;
  /* @} */

  /** \defgroup must This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*x*/, RangeType& /*ret*/) const = 0;
  /* @} */

  virtual RangeType evaluate(const DomainType& x) const;

  virtual void jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const;

  virtual JacobianRangeType jacobian(const DomainType& x) const;
}; // class FunctionInterface


// =====================================
// ===== LocalfunctionSetInterface =====
// =====================================
template <class E, class D, int d, class R, int r, int rC>
LocalfunctionSetInterface<E, D, d, R, r, rC>::LocalfunctionSetInterface(const EntityType& ent)
  : entity_(ent)
{
}

template <class E, class D, int d, class R, int r, int rC>
LocalfunctionSetInterface<E, D, d, R, r, rC>::~LocalfunctionSetInterface()
{
}

template <class E, class D, int d, class R, int r, int rC>
const typename LocalfunctionSetInterface<E, D, d, R, r, rC>::EntityType&
LocalfunctionSetInterface<E, D, d, R, r, rC>::entity() const
{
  return entity_;
}

template <class E, class D, int d, class R, int r, int rC>
std::vector<typename LocalfunctionSetInterface<E, D, d, R, r, rC>::RangeType>
LocalfunctionSetInterface<E, D, d, R, r, rC>::evaluate(const DomainType& xx) const
{
  std::vector<RangeType> ret(size(), RangeType(0));
  evaluate(xx, ret);
  return ret;
}

template <class E, class D, int d, class R, int r, int rC>
std::vector<typename LocalfunctionSetInterface<E, D, d, R, r, rC>::JacobianRangeType>
LocalfunctionSetInterface<E, D, d, R, r, rC>::jacobian(const DomainType& xx) const
{
  std::vector<JacobianRangeType> ret(size(), JacobianRangeType(0));
  jacobian(xx, ret);
  return ret;
}

template <class E, class D, int d, class R, int r, int rC>
bool LocalfunctionSetInterface<E, D, d, R, r, rC>::is_a_valid_point(const DomainType& xx) const
{
  const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity().type());
  return reference_element.checkInside(xx);
}


// ==================================
// ===== LocalfunctionInterface =====
// ==================================
template <class E, class D, int d, class R, int r, int rC>
LocalfunctionInterface<E, D, d, R, r, rC>::LocalfunctionInterface(const EntityType& ent)
  : BaseType(ent)
{
}

template <class E, class D, int d, class R, int r, int rC>
LocalfunctionInterface<E, D, d, R, r, rC>::~LocalfunctionInterface()
{
}

template <class E, class D, int d, class R, int r, int rC>
size_t LocalfunctionInterface<E, D, d, R, r, rC>::size() const
{
  return 1;
}

template <class E, class D, int d, class R, int r, int rC>
void LocalfunctionInterface<E, D, d, R, r, rC>::evaluate(const DomainType& xx, std::vector<RangeType>& ret) const
{
  assert(ret.size() >= 1);
  evaluate(xx, ret[0]);
}

template <class E, class D, int d, class R, int r, int rC>
void LocalfunctionInterface<E, D, d, R, r, rC>::jacobian(const DomainType& xx,
                                                         std::vector<JacobianRangeType>& ret) const
{
  assert(ret.size() >= 1);
  jacobian(xx, ret[0]);
}

template <class E, class D, int d, class R, int r, int rC>
typename LocalfunctionInterface<E, D, d, R, r, rC>::RangeType
LocalfunctionInterface<E, D, d, R, r, rC>::evaluate(const DomainType& xx) const
{
  RangeType ret(0);
  evaluate(xx, ret);
  return ret;
}

template <class E, class D, int d, class R, int r, int rC>
typename LocalfunctionInterface<E, D, d, R, r, rC>::JacobianRangeType
LocalfunctionInterface<E, D, d, R, r, rC>::jacobian(const DomainType& xx) const
{
  JacobianRangeType ret(0);
  jacobian(xx, ret);
  return ret;
}


// ========================================
// ===== LocalizableFunctionInterface =====
// ========================================
template <class E, class D, int d, class R, int r, int rC>
LocalizableFunctionInterface<E, D, d, R, r, rC>::~LocalizableFunctionInterface()
{
}

template <class E, class D, int d, class R, int r, int rC>
std::string LocalizableFunctionInterface<E, D, d, R, r, rC>::static_id()
{
  return "dune.stuff.function";
}

template <class E, class D, int d, class R, int r, int rC>
std::string LocalizableFunctionInterface<E, D, d, R, r, rC>::name() const
{
  return "dune.stuff.function";
}


// =============================
// ===== FunctionInterface =====
// =============================
template <class D, int d, class R, int r>
FunctionInterface<D, d, R, r>::~FunctionInterface()
{
}

template <class D, int d, class R, int r>
std::string FunctionInterface<D, d, R, r>::static_id()
{
  return "dune.stuff.function";
}

template <class D, int d, class R, int r>
std::string FunctionInterface<D, d, R, r>::name() const
{
  return "dune.stuff.function";
}

template <class D, int d, class R, int r>
int FunctionInterface<D, d, R, r>::order() const
{
  return -1;
}

template <class D, int d, class R, int r>
typename FunctionInterface<D, d, R, r>::RangeType FunctionInterface<D, d, R, r>::evaluate(const DomainType& x) const
{
  RangeType ret;
  evaluate(x, ret);
  return ret;
}

template <class D, int d, class R, int r>
void FunctionInterface<D, d, R, r>::jacobian(const DomainType& /*x*/, JacobianRangeType& /*ret*/) const
{
  DUNE_THROW(Dune::NotImplemented, "You really have to implement this!");
}

template <class D, int d, class R, int r>
typename FunctionInterface<D, d, R, r>::JacobianRangeType
FunctionInterface<D, d, R, r>::jacobian(const DomainType& x) const
{
  JacobianRangeType ret;
  jacobian(x, ret);
  return ret;
}


} // namespace Stuff
} // namespace Dune

#include "default.hh"

#define DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(etype, ddim)                                                        \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGE(Dune::Stuff::LocalfunctionSetInterface, etype, ddim)                     \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGE(Dune::Stuff::LocalfunctionInterface, etype, ddim)                        \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGE(Dune::Stuff::LocalizableFunctionInterface, etype, ddim)

#define DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGE(cname, etype, ddim)                                                \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGECOLS(cname, etype, ddim, 1)                                               \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGECOLS(cname, etype, ddim, 2)                                               \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGECOLS(cname, etype, ddim, 3)

#define DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGECOLS(cname, etype, ddim, rdim)                                      \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 1)                                     \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 2)                                     \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTION_INTERFACE_LIST_DOMAINFIELDTYPES(cname, etype, ddim, rdim, rcdim)                           \
  DUNE_STUFF_FUNCTION_INTERFACE_LIST_RANGEFIELDTYPES(cname, etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTION_INTERFACE_LIST_RANGEFIELDTYPES(cname, etype, dftype, ddim, rdim, rcdim)                    \
  DUNE_STUFF_FUNCTION_INTERFACE_LAST_EXPANSION(cname, etype, dftype, ddim, double, rdim, rcdim)                        \
  DUNE_STUFF_FUNCTION_INTERFACE_LAST_EXPANSION(cname, etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTION_INTERFACE_LAST_EXPANSION(cname, etype, dftype, ddim, rftype, rdim, rcdim)                  \
  extern template class cname<etype, dftype, ddim, rftype, rdim, rcdim>;

#include <dune/stuff/grid/fakeentity.hh>

typedef Dune::Stuff::Grid::FakeEntity<1> DuneStuffFake1dEntityType;
typedef Dune::Stuff::Grid::FakeEntity<2> DuneStuffFake2dEntityType;
typedef Dune::Stuff::Grid::FakeEntity<3> DuneStuffFake3dEntityType;

DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneStuffFake1dEntityType, 1)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneStuffFake2dEntityType, 2)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneStuffFake3dEntityType, 3)

#ifdef HAVE_DUNE_GRID

#include <dune/grid/sgrid.hh>

typedef typename Dune::SGrid<1, 1>::template Codim<0>::Entity DuneSGrid1dEntityType;
typedef typename Dune::SGrid<2, 2>::template Codim<0>::Entity DuneSGrid2dEntityType;
typedef typename Dune::SGrid<3, 3>::template Codim<0>::Entity DuneSGrid3dEntityType;

DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneSGrid3dEntityType, 3)

#include <dune/grid/yaspgrid.hh>

typedef typename Dune::YaspGrid<1>::template Codim<0>::Entity DuneYaspGrid1dEntityType;
typedef typename Dune::YaspGrid<2>::template Codim<0>::Entity DuneYaspGrid2dEntityType;
typedef typename Dune::YaspGrid<3>::template Codim<0>::Entity DuneYaspGrid3dEntityType;

DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#ifdef ALUGRID_CONFORM
#define DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#define ALUGRID_CONFORM 1
#endif
#ifdef ENABLE_ALUGRID
#define DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#define ENABLE_ALUGRID 1
#endif

#include <dune/grid/alugrid.hh>

typedef typename Dune::ALUSimplexGrid<2, 2>::template Codim<0>::Entity DuneAluSimplexGrid2dEntityType;
typedef typename Dune::ALUSimplexGrid<3, 3>::template Codim<0>::Entity DuneAluSimplexGrid3dEntityType;
typedef typename Dune::ALUCubeGrid<3, 3>::template Codim<0>::Entity DuneAluCubeGrid3dEntityType;

DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES(DuneAluCubeGrid3dEntityType, 3)

#ifdef DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTION_INTERFACE_ALUGRID_CONFORM_WAS_DEFINED_BEFORE
#else
#undef ALUGRID_CONFORM
#endif
#ifdef DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#undef DUNE_STUFF_FUNCTION_INTERFACE_ENABLE_ALUGRID_WAS_DEFINED_BEFORE
#else
#undef ENABLE_ALUGRID
#endif

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTION_INTERFACE_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTION_INTERFACE_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTION_INTERFACE_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTION_INTERFACE_LIST_DIMRANGE
#undef DUNE_STUFF_FUNCTION_INTERFACE_LIST_CLASSES

extern template class Dune::Stuff::FunctionInterface<double, 1, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 1, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 1, double, 3>;

extern template class Dune::Stuff::FunctionInterface<double, 2, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 2, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 2, double, 3>;

extern template class Dune::Stuff::FunctionInterface<double, 3, double, 1>;
extern template class Dune::Stuff::FunctionInterface<double, 3, double, 2>;
extern template class Dune::Stuff::FunctionInterface<double, 3, double, 3>;

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
