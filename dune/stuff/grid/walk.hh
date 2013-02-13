#ifndef DUNE_STUFF_WALK_HH_INCLUDED
#define DUNE_STUFF_WALK_HH_INCLUDED

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/static_assert.hh>
#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/grid/common/geometry.hh>
#include <vector>
#include <boost/format.hpp>

namespace Dune {
namespace Stuff {
namespace Grid {

/** \brief Useful dummy functor if you don't have anything to do on entities/intersections
 **/
struct GridWalkDummyFunctor
{
  GridWalkDummyFunctor()
  {
  }

  template <class Entity>
  void operator()(const Entity&, const int) const
  {
  }
  template <class Entity, class Intersection>
  void operator()(const Entity&, const Intersection&) const
  {
  }
};

namespace {
/** \brief global \ref GridWalkDummyFunctor instance
 **/
const GridWalkDummyFunctor gridWalkDummyFunctor;
}

/** \brief applies Functors on each \ref Entity/\ref Intersection of a given \ref GridView
 *  \todo allow stacking of functor to save gridwalks?
 *  \tparam GridViewImp any \ref GridView interface compliant type
 *  \tparam codim determines the codim of the Entities that are iterated on
 **/
template <class GridViewImp, int codim = 0>
class GridWalk
{
  typedef Dune::GridView<typename GridViewImp::Traits> GridViewType;

public:
  GridWalk(const GridViewType& gp)
    : gridView_(gp)
  {
  }

  /** \param entityFunctor is applied on all codim 0 entities presented by \var gridView_
   *  \param intersectionFunctor is applied on all Intersections of all codim 0 entities
   *        presented by \var gridView_
   *  \note only instantiable for codim == 0
   */
  template <class EntityFunctor, class IntersectionFunctor>
  void operator()(EntityFunctor& entityFunctor, IntersectionFunctor& intersectionFunctor) const
  {
    dune_static_assert(codim == 0, "walking intersections is only possible for codim 0 entities");
    for (const auto& entity : Common::viewRange(gridView_)) {
      const int entityIndex = gridView_.indexSet().index(entity);
      entityFunctor(entity, entityIndex);
      for (const auto& intersection : DSC::intersectionRange(gridView_, entity)) {
        intersectionFunctor(entity, intersection);
      }
    }
  }

  /** \param entityFunctor is applied on all codim entities presented by \var gridView_
   *  \note only instantiable for codim < GridView::dimension
   */
  template <class EntityFunctor>
  void operator()(EntityFunctor& entityFunctor) const
  {
    dune_static_assert(codim <= GridViewType::dimension, "codim too high to walk");
    for (const auto& entity : DSC::viewRange(gridView_)) {
      const int entityIndex = gridView_.indexSet().index(entity);
      entityFunctor(entity, entityIndex);
    }
  }

  template <class Functor>
  void walkCodim0(Functor& f) const DUNE_DEPRECATED_MSG("use operator()(Functor) instead ");


private:
  const GridViewType& gridView_;
};

template <class V, int i>
template <class Functor>
void GridWalk<V, i>::walkCodim0(Functor& f) const
{
  this->operator()(f);
}

//!
template <class ViewImp, int codim = 0>
GridWalk<Dune::GridView<ViewImp>, codim> make_gridwalk(const Dune::GridView<ViewImp>& view)
{
  return GridWalk<Dune::GridView<ViewImp>, codim>(view);
}

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // ifndef DUNE_STUFF_WALK_HH_INCLUDED
/** Copyright (c) 2012, Felix Albrecht, Rene Milk
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/
