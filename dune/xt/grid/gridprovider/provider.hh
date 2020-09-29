// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2019)
//   Kirsten Weber   (2012)
//   René Fritze     (2012 - 2019)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016 - 2017, 2020)

#ifndef DUNE_XT_GRID_PROVIDER_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_PROVIDER_HH

#include <memory>
#include <set>
#include <string>

#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/capabilities.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/output/entity_visualization.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class GridImp>
class GridProvider
{
  static_assert(is_grid<GridImp>::value, "");

public:
  using ThisType = GridProvider<GridImp>;

  using GridType = GridImp;
  static constexpr size_t dimDomain = GridImp::dimension;
  using DomainFieldType = typename GridType::ctype;
  using DomainType = FieldVector<DomainFieldType, dimDomain>;
  using EntityType = typename GridType::template Codim<0>::Entity;
  using LevelGridViewType = typename Layer<GridType, Layers::level, Backends::view>::type;
  using LeafGridViewType = typename Layer<GridType, Layers::leaf, Backends::view>::type;


  static const std::string static_id()
  {
    return "xt.grid.gridprovider";
  }

  /**
   * \attention Do not delete grd_ptr manually afterwards!
   */
  GridProvider(GridType*&& grd_ptr)
    : grid_ptr_(std::move(grd_ptr))
  {}

  GridProvider(std::shared_ptr<GridType> grd_ptr)
    : grid_ptr_(grd_ptr)
  {}

#if DUNE_VERSION_GTE(DUNE_COMMON, 2, 7)
  GridProvider(Dune::ToUniquePtr<GridType> grd_ptr)
    : grid_ptr_(std::shared_ptr<GridType>(grd_ptr))
  {}
#endif

  GridProvider(const ThisType& other) = default;

  // Manual ctor required for clang 3.8.1-12~bpo8+1 (otherwise: undefined reference).
  GridProvider(ThisType&& source)
    : grid_ptr_(source.grid_ptr_)
  {}

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const std::shared_ptr<GridType>& grid_ptr() const
  {
    return grid_ptr_;
  }

  std::shared_ptr<GridType> grid_ptr()
  {
    return grid_ptr_;
  }

  const GridType& grid() const
  {
    return *grid_ptr_;
  }

  GridType& grid()
  {
    return *grid_ptr_;
  }

  int max_level() const
  {
    return grid_ptr_->maxLevel();
  }

  template <Layers lr, Backends backend>
  typename Layer<GridType, lr, backend>::type layer(const int lvl = 0) const
  {
    return Layer<GridType, lr, backend>::create(*grid_ptr_, lvl);
  }

  LevelGridViewType level_view(const int lvl) const
  {
    return Layer<GridType, Layers::level, Backends::view>::create(*grid_ptr_, lvl);
  }

  LeafGridViewType leaf_view() const
  {
    return Layer<GridType, Layers::leaf, Backends::view>::create(*grid_ptr_);
  }

  void visualize(const std::string filename = static_id(),
                 const Common::Configuration& boundary_info_cfg = Common::Configuration()) const
  {
    if (boundary_info_cfg.empty())
      visualize_plain(filename);
    else
      visualize_with_boundary(boundary_info_cfg, filename);
  } // ... visualize(...)

  void visualize(const Common::Configuration& boundary_info_cfg) const
  {
    visualize_with_boundary(boundary_info_cfg, static_id());
  }

  void visualize(const Common::Configuration& boundary_info_cfg, const std::string filename) const
  {
    visualize_with_boundary(boundary_info_cfg, filename);
  }

private:
  template <class G>
  struct global_refine_helper
  {
    void operator()(G& g, int count)
    {
      g.preAdapt();
      g.globalRefine(count);
      g.postAdapt();
      g.loadBalance();
    }
  }; // struct refine_helper

#if HAVE_ALBERTA

  template <int d, int dW>
  struct global_refine_helper<AlbertaGrid<d, dW>>
  {
    using G = AlbertaGrid<d, dW>;

    void operator()(G& g, int count)
    {
      g.globalRefine(count);
      g.loadBalance();
    }
  }; // struct refine_helper

#endif // HAVE_ALBERTA

public:
  void global_refine(const int count = 1)
  {
    global_refine_helper<GridType>()(grid(), count);
  }

private:
  void visualize_plain(const std::string filename) const
  {
    if (GridType::dimension > 3) // give us a call if you have any idea!
      DUNE_THROW(NotImplemented, "For grids of dimension > 3!");
    for (auto lvl : Common::value_range(max_level() + 1)) {
      auto grid_view = level_view(lvl);
      // vtk writer
      Dune::VTKWriter<LevelGridViewType> vtkwriter(grid_view);
      // codim 0 entity id
      const auto entityId = ElementVisualization::IndexFunctor<LevelGridViewType>(grid_view).values(grid_view);
      vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
      // boundary id
      const std::vector<double> boundary_ids =
          ElementVisualization::BoundaryIDFunctor<LevelGridViewType>(grid_view).values(grid_view);
      vtkwriter.addCellData(boundary_ids, "boundary_id__level_" + Common::to_string(lvl));
      // write
      vtkwriter.write(filename + "__level_" + Common::to_string(lvl), VTK::appendedraw);
    }
  } // ... visualize_plain(...)

  void visualize_with_boundary(const Common::Configuration& boundary_info_cfg, const std::string filename) const
  {
    if (GridType::dimension > 3) // give us a call if you have any idea!
      DUNE_THROW(NotImplemented, "For grids of dimension > 3!");
    // boundary info
    using BoundaryInfoFactory = XT::Grid::BoundaryInfoFactory<typename LevelGridViewType::Intersection>;
    auto boundary_info_ptr = BoundaryInfoFactory::create(boundary_info_cfg);
    for (auto lvl : Common::value_range(max_level() + 1)) {
      auto grid_view = level_view(lvl);
      // vtk writer
      Dune::VTKWriter<LevelGridViewType> vtkwriter(grid_view);
      // codim 0 entity id
      const auto entityId = ElementVisualization::IndexFunctor<LevelGridViewType>(grid_view).values(grid_view);
      vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
      // boundary id
      const auto boundary_ids = ElementVisualization::BoundaryIDFunctor<LevelGridViewType>(grid_view).values(grid_view);
      vtkwriter.addCellData(boundary_ids, "boundary_id__level_" + Common::to_string(lvl));
      // dirichlet values
      const auto dirichlet = ElementVisualization::BoundaryTypeFunctor<LevelGridViewType, decltype(*boundary_info_ptr)>(
                                 grid_view, *boundary_info_ptr, "dirichlet")
                                 .values(grid_view);
      vtkwriter.addCellData(dirichlet, "isDirichletBoundary__level_" + Common::to_string(lvl));
      // neumann values
      const auto neumann = ElementVisualization::BoundaryTypeFunctor<LevelGridViewType, decltype(*boundary_info_ptr)>(
                               grid_view, *boundary_info_ptr, "neumann")
                               .values(grid_view);
      vtkwriter.addCellData(neumann, "isNeumannBoundary__level_" + Common::to_string(lvl));
      // write
      vtkwriter.write(filename + "__level_" + Common::to_string(lvl), VTK::appendedraw);
    }
  } // ... visualize_with_boundary(...)

  std::shared_ptr<GridType> grid_ptr_;
}; // class GridProvider


} // namespace Grid
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_GRID_PROVIDER_PROVIDER_HH
