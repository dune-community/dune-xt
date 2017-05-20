// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016)

#ifndef DUNE_XT_GRID_PROVIDER_EOC_HH
#define DUNE_XT_GRID_PROVIDER_EOC_HH

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/dgfparser.hh>
#endif

#include <dune/xt/grid/gridprovider/provider.hh>

namespace Dune {
namespace XT {
namespace Grid {

/**
 *  The purpose of this class is to behave like a XT::Grid::GridProvider and at the same time to provide a
 *  means to obtain the real grid level corresponding to a refinement level.
 */
template <class GridImp>
class EOCGridProvider : public GridProvider<GridImp>
{
  typedef GridProvider<GridImp> BaseType;

public:
  using typename BaseType::GridType;

  EOCGridProvider(const Common::Configuration& grid_cfg, const size_t num_refs)
    : BaseType(GridProviderFactory<GridType>::create(grid_cfg).grid_ptr())
  {
    setup(num_refs);
  }

  size_t num_refinements() const
  {
    assert(levels_.size() > 0);
    return levels_.size() - 1;
  }

  int level_of(const size_t refinement) const
  {
    assert(refinement <= num_refinements());
    return levels_[refinement];
  }

  int reference_level() const
  {
    return reference_level_;
  }

  typename Layer<GridType, Layers::level, Backends::view>::type reference_grid_view() const
  {
    return this->level_view(reference_level_);
  }

private:
  void setup(const size_t num_refinements)
  {
    levels_.push_back(this->grid().maxLevel());
    static const int refine_steps_for_half = DGFGridInfo<GridType>::refineStepsForHalf();
    for (size_t rr = 0; rr < num_refinements; ++rr) {
      this->grid().globalRefine(refine_steps_for_half);
      levels_.push_back(this->grid().maxLevel());
    }
    this->grid().globalRefine(refine_steps_for_half);
    reference_level_ = this->grid().maxLevel();
  } // ... setup(...)

  std::vector<int> levels_;
  int reference_level_;
}; // class EOCGridProvider


} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_EOC_HH
