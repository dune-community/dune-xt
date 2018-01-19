// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016, 2018)

#ifndef DUNE_XT_GRID_PROVIDER_EOC_HH
#define DUNE_XT_GRID_PROVIDER_EOC_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/gridprovider/factory.hh>

namespace Dune {
namespace XT {
namespace Grid {


/**
 *  The purpose of this class is to behave like a XT::Grid::GridProvider and at the same time to provide a
 *  means to obtain the real grid level corresponding to a refinement level.
 */
template <class GridImp>
class LevelBasedEOCGridProvider : public GridProvider<GridImp>
{
  typedef GridProvider<GridImp> BaseType;

public:
  using typename BaseType::GridType;

  static const constexpr Layers layer_type = Layers::level;

  LevelBasedEOCGridProvider(const Common::Configuration& grid_cfg, const size_t num_refs)
    : BaseType(GridProviderFactory<GridType>::create(grid_cfg).grid_ptr())
  {
    levels_.push_back(this->grid().maxLevel());
    const auto refine_steps_for_half = DGFGridInfo<GridType>::refineStepsForHalf();
    for (size_t rr = 0; rr < num_refs; ++rr) {
      this->grid().globalRefine(refine_steps_for_half);
      levels_.push_back(this->grid().maxLevel());
    }
    this->grid().globalRefine(refine_steps_for_half);
    reference_level_ = this->grid().maxLevel();
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

  const BaseType& level_provider(const size_t /*refinement*/) const
  {
    return *this;
  }

  BaseType& level_provider(const size_t /*refinement*/)
  {
    return *this;
  }

  const BaseType& reference_provider() const
  {
    return *this;
  }

  BaseType& reference_provider()
  {
    return *this;
  }

  typename Layer<GridType, Layers::level, Backends::view>::type reference_grid_view() const
  {
    return this->level_view(reference_level_);
  }

private:
  std::vector<int> levels_;
  int reference_level_;
}; // class LevelBasedEOCGridProvider


template <class GridImp>
class LeafBasedEOCGridProvider
{
  typedef GridProvider<GridImp> GridProviderType;

public:
  typedef typename GridProviderType::GridType GridType;
  typedef typename GridProviderType::LeafGridViewType LevelGridViewType;

  static const constexpr Layers layer_type = Layers::leaf;

  LeafBasedEOCGridProvider(const Common::Configuration& grid_cfg, const size_t num_refs)
  {
    level_grids_.emplace_back(new GridProviderType(GridProviderFactory<GridType>::create(grid_cfg)));
    const auto refine_steps_for_half = DGFGridInfo<GridType>::refineStepsForHalf();
    for (size_t rr = 0; rr < num_refs; ++rr) {
      level_grids_.emplace_back(new GridProviderType(GridProviderFactory<GridType>::create(grid_cfg)));
      level_grids_.back()->grid().globalRefine((rr + 1) * refine_steps_for_half);
    }
    reference_grid_ = std::make_unique<GridProviderType>(GridProviderFactory<GridType>::create(grid_cfg));
    reference_grid_->grid().globalRefine((num_refs + 1) * refine_steps_for_half);
  }

  size_t num_refinements() const
  {
    assert(level_grids_.size() > 0);
    return level_grids_.size() - 1;
  }

  int level_of(const size_t refinement) const
  {
    return -1;
  }

  int reference_level() const
  {
    return -1;
  }

  const GridProviderType& level_provider(const size_t refinement) const
  {
    return *level_grids_.at(refinement);
  }

  GridProviderType& level_provider(const size_t refinement)
  {
    return *level_grids_.at(refinement);
  }

  const GridProviderType& reference_provider() const
  {
    return *reference_grid_;
  }

  GridProviderType& reference_provider()
  {
    return *reference_grid_;
  }

  LevelGridViewType reference_grid_view() const
  {
    return reference_grid_->leaf_view();
  }

private:
  std::vector<std::unique_ptr<GridProviderType>> level_grids_;
  std::unique_ptr<GridProviderType> reference_grid_;
}; // class LeafBasedEOCGridProvider


template <class GridImp>
class DdSubdomainsBasedEOCGridProvider
{
  typedef GridProvider<GridImp, DD::SubdomainGrid<GridImp>> GridProviderType;
  typedef DdSubdomainGridProviderFactory<GridImp> GridProviderFactoryType;

public:
  typedef typename GridProviderType::GridType GridType;
  typedef typename GridProviderType::LeafGridViewType LevelGridViewType;

  static const constexpr Layers layer_type = Layers::leaf;

  DdSubdomainsBasedEOCGridProvider(const Common::Configuration& grid_cfg, const size_t num_refs)
  {
    const auto refine_steps_for_half = DGFGridInfo<GridType>::refineStepsForHalf();
    Common::Configuration cfg = grid_cfg;
    level_grids_.emplace_back(new GridProviderType(GridProviderFactoryType::create(cfg)));
    for (size_t rr = 0; rr < num_refs; ++rr) {
      cfg["num_refinements"] = XT::Common::to_string(cfg.get<size_t>("num_refinements") + refine_steps_for_half);
      level_grids_.emplace_back(new GridProviderType(GridProviderFactoryType::create(cfg)));
    }
    cfg["num_refinements"] = XT::Common::to_string(cfg.get<size_t>("num_refinements") + refine_steps_for_half);
    reference_grid_ = std::make_unique<GridProviderType>(GridProviderFactoryType::create(cfg));
  }

  size_t num_refinements() const
  {
    assert(level_grids_.size() > 0);
    return level_grids_.size() - 1;
  }

  int level_of(const size_t refinement) const
  {
    return -1;
  }

  int reference_level() const
  {
    return -1;
  }

  const GridProviderType& level_provider(const size_t refinement) const
  {
    return *level_grids_.at(refinement);
  }

  GridProviderType& level_provider(const size_t refinement)
  {
    return *level_grids_.at(refinement);
  }

  const GridProviderType& reference_provider() const
  {
    return *reference_grid_;
  }

  GridProviderType& reference_provider()
  {
    return *reference_grid_;
  }

  LevelGridViewType reference_grid_view() const
  {
    return reference_grid_->leaf_view();
  }

private:
  std::vector<std::unique_ptr<GridProviderType>> level_grids_;
  std::unique_ptr<GridProviderType> reference_grid_;
}; // class DdSubdomainsBasedEOCGridProvider


template <class G, class DdGrid = int>
class EOCGridProvider : public LevelBasedEOCGridProvider<G>
{
public:
  template <class... Args>
  EOCGridProvider(Args&&... args)
    : LevelBasedEOCGridProvider<G>(std::forward<Args>(args)...)
  {
  }
};

template <class G>
class EOCGridProvider<G, DD::SubdomainGrid<G>> : public DdSubdomainsBasedEOCGridProvider<G>
{
public:
  template <class... Args>
  EOCGridProvider(Args&&... args)
    : DdSubdomainsBasedEOCGridProvider<G>(std::forward<Args>(args)...)
  {
  }
};

#if HAVE_DUNE_ALUGRID

template <class Comm>
class EOCGridProvider<Dune::ALUGrid<2, 2, simplex, conforming, Comm>>
    : public LeafBasedEOCGridProvider<Dune::ALUGrid<2, 2, simplex, conforming, Comm>>
{
public:
  template <class... Args>
  EOCGridProvider(Args&&... args)
    : LeafBasedEOCGridProvider<Dune::ALUGrid<2, 2, simplex, conforming, Comm>>(std::forward<Args>(args)...)
  {
  }
};

#endif // HAVE_DUNE_ALUGRID


} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_EOC_HH
