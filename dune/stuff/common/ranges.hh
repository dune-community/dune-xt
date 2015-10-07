// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_RANGES_RANGES_HH
#define DUNE_STUFF_RANGES_RANGES_HH

#include <type_traits>

#include <boost/serialization/static_warning.hpp>
#include <boost/iterator/iterator_facade.hpp>

#if HAVE_DUNE_GRID
#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entity.hh>
#endif // HAVE_DUNE_GRID

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>
#endif // HAVE_DUNE_FEM

#include <dune/stuff/common/math.hh>

namespace Dune {

#if HAVE_DUNE_FEM

namespace Fem {

template <class DiscreteFunctionTraits>
auto begin(const Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& func) -> decltype(func.dbegin())
{
  return func.dbegin();
}

template <class DiscreteFunctionTraits>
auto end(const Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& func) -> decltype(func.dend())
{
  return func.dend();
}

template <class DiscreteFunctionTraits>
auto begin(Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& func) -> decltype(func.dbegin())
{
  return func.dbegin();
}

template <class DiscreteFunctionTraits>
auto end(Dune::Fem::DiscreteFunctionInterface<DiscreteFunctionTraits>& func) -> decltype(func.dend())
{
  return func.dend();
}

} // namespace Fem

#endif // HAVE_DUNE_FEM

namespace Stuff {
namespace Common {

#if HAVE_DUNE_GRID

//! adapter enabling view usage in range-based for
template <class GridPartViewType, size_t codim = 0>
class EntityRange
{
  const GridPartViewType& part_view_;

public:
  EntityRange(const GridPartViewType& part_view)
    : part_view_(part_view)
  {
  }

  auto begin() const -> decltype(part_view_.template begin<codim>())
  {
    return part_view_.template begin<codim>();
  }

  auto end() const -> decltype(part_view_.template end<codim>())
  {
    return part_view_.template end<codim>();
  }
}; // class EntityRange

template <class GridViewTraits, size_t codim = 0>
EntityRange<Dune::GridView<GridViewTraits>, codim> entityRange(const Dune::GridView<GridViewTraits>& view)
{
  return EntityRange<Dune::GridView<GridViewTraits>, codim>(view);
}

#if HAVE_DUNE_FEM

template <class GP, size_t codim = 0>
EntityRange<Dune::Fem::GridPartInterface<GP>, codim> entityRange(const Dune::Fem::GridPartInterface<GP>& part)
{
  return EntityRange<Dune::Fem::GridPartInterface<GP>, codim>(part);
}

#endif // HAVE_DUNE_FEM

/** adapter enabling intersectionniterator usage in range-based for
 * works for GridParts and GridViews
 */
template <class GridAbstractionType, class EntityType>
class IntersectionRange
{
  const GridAbstractionType& view_;
  const EntityType& entity_;

public:
  IntersectionRange(const GridAbstractionType& view, const EntityType& entity)
    : view_(view)
    , entity_(entity)
  {
  }

  auto begin() const -> decltype(view_.ibegin(entity_))
  {
    return view_.ibegin(entity_);
  }

  auto end() const -> decltype(view_.iend(entity_))
  {
    return view_.iend(entity_);
  }
};

template <class GridViewTraits>
IntersectionRange<Dune::GridView<GridViewTraits>, typename Dune::GridView<GridViewTraits>::template Codim<0>::Entity>
intersectionRange(const Dune::GridView<GridViewTraits>& gridview,
                  const typename Dune::GridView<GridViewTraits>::template Codim<0>::Entity& entity)
{
  return IntersectionRange<Dune::GridView<GridViewTraits>,
                           typename Dune::GridView<GridViewTraits>::template Codim<0>::Entity>(gridview, entity);
}

//! custom const iterator for \ref FixedMap
template <class GeometryType>
class ConstCornerIterator
    : public boost::iterator_facade<ConstCornerIterator<GeometryType>, const typename GeometryType::GlobalCoordinate,
                                    boost::forward_traversal_tag>
{
  typedef ConstCornerIterator<GeometryType> ThisType;
  typedef typename GeometryType::GlobalCoordinate CornerType;

public:
  ConstCornerIterator()
    : index_(-1)
    , geometry_(nullptr)
  {
  }

  explicit ConstCornerIterator(const GeometryType* const geometry, int i = 0)
    : index_(i)
    , geometry_(geometry)
  {
  }

private:
  friend class boost::iterator_core_access;

  void increment()
  {
    index_++;
  }

  bool equal(ThisType const& other) const
  {
    return this->geometry_ && (index_ == other.index_) && (this->geometry_ == other.geometry_);
  }

  CornerType dereference() const
  {
    return geometry_->corner(index_);
  }

  int index_;
  const GeometryType* const geometry_;
};

template <class GeometryType>
class CornerRange
{
  typedef ConstCornerIterator<GeometryType> IteratorType;

public:
  CornerRange(const GeometryType& geometry)
    : geometry_(geometry)
  {
  }

  IteratorType begin() const
  {
    return IteratorType(&geometry_, 0);
  }

  IteratorType end() const
  {
    return IteratorType(&geometry_, geometry_.corners());
  }

private:
  const GeometryType& geometry_;
};

template <int mydim, int cdim, class GridImp, template <int, int, class> class GeometryImp>
CornerRange<Dune::Geometry<mydim, cdim, GridImp, GeometryImp>>
cornerRange(const Dune::Geometry<mydim, cdim, GridImp, GeometryImp>& geometry)
{
  return CornerRange<Dune::Geometry<mydim, cdim, GridImp, GeometryImp>>(geometry);
}

template <int mydim, int cdim, class GridImp, template <int, int, class> class EntityImp>
auto cornerRange(const Dune::Entity<mydim, cdim, GridImp, EntityImp>& entity)
    -> CornerRange<typename std::remove_reference<decltype(entity.geometry())>::type>
{
  return CornerRange<typename std::remove_reference<decltype(entity.geometry())>::type>(entity.geometry());
}

#endif //#if HAVE_DUNE_GRID

#if HAVE_DUNE_FEM

//! Range adapter for lagrange points from lagrange spaces
template <class GridPartType, int order, int faceCodim>
class LagrangePointSetRange
{
  typedef Dune::Fem::LagrangePointSet<GridPartType, order> LagrangePointSetType;
  typedef typename LagrangePointSetType::template Codim<faceCodim>::SubEntityIteratorType SubEntityIteratorType;
  const LagrangePointSetType& lp_set_;
  const size_t subEntity_;

public:
  /** the template isn't lazyness here, the underlying set is templated on it too
   */
  template <class DiscreteFunctionspaceType, class EntityType>
  LagrangePointSetRange(const DiscreteFunctionspaceType& space, const EntityType& entity, const size_t subEntity)
    : lp_set_(space.lagrangePointSet(entity))
    , subEntity_(subEntity)
  {
  }

  LagrangePointSetRange(const LagrangePointSetType& lp_set, const size_t subEntity)
    : lp_set_(lp_set)
    , subEntity_(subEntity)
  {
  }

  SubEntityIteratorType begin() const
  {
    return lp_set_.template beginSubEntity<faceCodim>(subEntity_);
  }
  SubEntityIteratorType end() const
  {
    return lp_set_.template endSubEntity<faceCodim>(subEntity_);
  }
};

template <size_t codim, class DiscreteFunctionspaceType, class EntityType>
LagrangePointSetRange<typename DiscreteFunctionspaceType::GridPartType, DiscreteFunctionspaceType::polynomialOrder,
                      codim>
lagrangePointSetRange(const DiscreteFunctionspaceType& space, const EntityType& entity, const size_t subEntity)
{
  return LagrangePointSetRange<typename DiscreteFunctionspaceType::GridPartType,
                               DiscreteFunctionspaceType::polynomialOrder,
                               codim>(space, entity, subEntity);
}

template <class LgPointSetType, size_t codim = 1>
LagrangePointSetRange<typename LgPointSetType::GridPartType, LgPointSetType::polynomialOrder, codim>
lagrangePointSetRange(const LgPointSetType& lpset, const size_t subEntity)
{
  return LagrangePointSetRange<typename LgPointSetType::GridPartType, LgPointSetType::polynomialOrder, codim>(
      lpset, subEntity);
}

template <class GridPartTraits>
IntersectionRange<Dune::Fem::GridPartInterface<GridPartTraits>,
                  typename Dune::Fem::GridPartInterface<GridPartTraits>::template Codim<0>::EntityType>
intersectionRange(const Dune::Fem::GridPartInterface<GridPartTraits>& gridpart,
                  const typename Dune::Fem::GridPartInterface<GridPartTraits>::template Codim<0>::EntityType& entity)
{
  return IntersectionRange<Dune::Fem::GridPartInterface<GridPartTraits>,
                           typename Dune::Fem::GridPartInterface<GridPartTraits>::template Codim<0>::EntityType>(
      gridpart, entity);
}

#endif // HAVE_DUNE_FEM

//! get a vector with values in [start : increment : end)
template <class T, class sequence = std::vector<T>>
typename std::enable_if<!std::is_enum<T>::value, sequence>::type valueRange(const T start, const T end,
                                                                            const T increment = Epsilon<T>::value)
{
  // sadly, no overloaded version of std::abs is available for
  // unsigned long long, so we compute the absolute value of increment
  // ourselves
  const auto incrementAbs = Dune::Stuff::Common::abs(increment);
  assert(incrementAbs > 0);
  sequence ret(typename sequence::size_type(((end > start) ? end - start : start - end) / incrementAbs), start);
  typename sequence::size_type i = 0;
  std::generate(std::begin(ret), std::end(ret), [&]() { return T(start + (increment * i++)); });
  return ret;
}

//! signature for enumeration Types T
template <class T, class sequence = std::vector<typename absretval<T>::type>>
typename std::enable_if<std::is_enum<T>::value, sequence>::type
valueRange(const T start, const T end,
           const typename absretval<T>::type increment = Epsilon<typename absretval<T>::type>::value)
{
  typedef typename absretval<T>::type R;
  return valueRange(static_cast<R>(start), static_cast<R>(end), increment);
}

//! get a vector with values in [0 : Epsilon<T> : end)
template <class T, class sequence = std::vector<T>>
typename std::enable_if<!std::is_enum<T>::value, sequence>::type valueRange(const T end)
{
  return valueRange(T(0), end);
}

//! get a vector with values in [0 : Epsilon<T> : end)
template <class T, class sequence = std::vector<typename absretval<T>::type>>
typename std::enable_if<std::is_enum<T>::value, sequence>::type valueRange(const T end)
{
  return valueRange(T(0), end);
}

} // namespace Common
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_RANGES_RANGES_HH
