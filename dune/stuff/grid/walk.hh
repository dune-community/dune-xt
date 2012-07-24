#ifndef DUNE_STUFF_WALK_HH_INCLUDED
#define DUNE_STUFF_WALK_HH_INCLUDED

#include <dune/common/static_assert.hh>
#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>
#include <dune/stuff/common/math.hh>
#include <dune/stuff/common/misc.hh>
#include <dune/grid/common/geometry.hh>
#include <vector>
#include <boost/format.hpp>

namespace Dune {
namespace Stuff {
namespace Grid {

// ! Base class for Gridwalk Functors that don't want to reimplement pre/postWalk
struct WalkFunctorDefault
{
  void preWalk() const
  {
  }
  void postWalk() const
  {
  }
};

/** \brief lets you apply a Functor to each entity
   * \todo allow stacking of operators to save gridwalks
   * \todo threadsafe maps (haha:P)
   */
template <class GridView, int codim = 0>
class Walk
{
private:
  typedef typename GridView::template Codim<0>::Iterator ElementIterator;
  typedef typename GridView::IntersectionIterator IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection::EntityPointer EntityPointer;

public:
  Walk(const GridView& gp)
    : gridView_(gp)
  {
  }

  template <class Functor>
  void operator()(Functor& f) const
  {
    f.preWalk();
    for (ElementIterator it = gridView_.template begin<0>(); it != gridView_.template end<0>(); ++it) {
      const int ent_idx = gridView_.indexSet().index(*it);
      f(*it, *it, ent_idx, ent_idx);
      IntersectionIteratorType intItEnd = gridView_.iend(*it);
      for (IntersectionIteratorType intIt = gridView_.ibegin(*it); intIt != intItEnd; ++intIt) {
        const auto& intersection = *intIt;
        if (!intersection.boundary()) {
          const auto neighbour_ptr = intersection.outside();
          const int neigh_idx = gridView_.indexSet().index(*neighbour_ptr);
          f(*it, *neighbour_ptr, ent_idx, neigh_idx);
        }
      }
    }
    f.postWalk();
  } // ()

  template <class Functor>
  void walkCodim0(Functor& f) const
  {
    for (ElementIterator it = gridView_.template begin<0>(); it != gridView_.template end<0>(); ++it) {
      const int ent_idx = gridView_.indexSet().index(*it);
      f(*it, ent_idx);
    }
  } // walkCodim0

private:
  const GridView& gridView_;
};

} // namespace Grid
} // namespace Stud
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
