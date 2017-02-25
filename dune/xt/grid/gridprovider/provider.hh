// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2012)
//   Rene Milk       (2012 - 2016)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_XT_GRID_PROVIDER_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_PROVIDER_HH

#include <memory>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/ranges.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/rangegenerators.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/layers.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class GridImp>
class GridProvider
{
  static_assert(is_grid<GridImp>::value, "");

public:
  typedef GridProvider<GridImp> ThisType;

  typedef GridImp GridType;
  static const size_t dimDomain = GridImp::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef typename Layer<GridType, Layers::level, Backends::view>::type LevelGridViewType;
#if HAVE_DUNE_FEM
  typedef typename Layer<GridType, Layers::level, Backends::part>::type LevelGridPartType;
#endif

  typedef typename Layer<GridType, Layers::leaf, Backends::view>::type LeafGridViewType;
#if HAVE_DUNE_FEM
  typedef typename Layer<GridType, Layers::leaf, Backends::part>::type LeafGridPartType;
  typedef typename Layer<GridType, Layers::adaptive_leaf, Backends::part>::type AdaptiveLeafGridPartType;
#endif

  static const std::string static_id()
  {
    return "xt.grid.gridprovider";
  }

  /**
   * \attention Do not delete grd_ptr manually afterwards!
   */
  GridProvider(GridType* grd_ptr)
    : grid_ptr_(grd_ptr)
  {
  }

  GridProvider(std::unique_ptr<GridType>&& grd_ptr)
    : grid_ptr_(grd_ptr)
  {
  }

  GridProvider(std::shared_ptr<GridType> grd_ptr)
    : grid_ptr_(grd_ptr)
  {
  }

  GridProvider(const ThisType& other) = default;
  GridProvider(ThisType&& source) = default;

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

  template <Backends backend>
  typename Layer<GridType, Layers::level, backend>::type level(const int lvl = 0) const
  {
    return Layer<GridType, Layers::level, backend>::create(*grid_ptr_, lvl);
  }

  template <Backends backend>
  typename Layer<GridType, Layers::leaf, backend>::type leaf() const
  {
    return Layer<GridType, Layers::leaf, backend>::create(*grid_ptr_);
  }

  template <Layers lr, Backends backend>
  typename Layer<GridType, lr, backend>::type layer(const int lvl = 0) const
  {
    return Layer<GridType, lr, backend>::create(*grid_ptr_, lvl);
  }

  template <Layers lr>
  typename Layer<GridType, lr, Backends::part>::type layer(const int lvl = 0)
  {
    return Layer<GridType, lr, Backends::part>::create(*grid_ptr_, lvl);
  }

  LevelGridViewType level_view(const int lvl) const
  {
    return Layer<GridType, Layers::level, Backends::view>::create(*grid_ptr_, lvl);
  }

  LeafGridViewType leaf_view() const
  {
    return Layer<GridType, Layers::leaf, Backends::view>::create(*grid_ptr_);
  }

#if HAVE_DUNE_FEM

  LevelGridPartType level_part(const int lvl)
  {
    return Layer<GridType, Layers::level, Backends::part>::create(*grid_ptr_, lvl);
  }

  LeafGridPartType leaf_part()
  {
    return Layer<GridType, Layers::leaf, Backends::part>::create(*grid_ptr_);
  }

  AdaptiveLeafGridPartType adaptive_leaf_part()
  {
    return Layer<GridType, Layers::adaptive_leaf, Backends::part>::create(*grid_ptr_);
  }

#endif // HAVE_DUNE_FEM

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
  template <class G, bool anything = true>
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

  template <int d, int dW, bool anything>
  struct global_refine_helper<AlbertaGrid<d, dW>, anything>
  {
    typedef AlbertaGrid<d, dW> G;

    void operator()(G& g, int count)
    {
      g.globalRefine(count);
      g.loadBalance();
    }
  }; // struct refine_helper

#endif // HAVE_ALBERTA

public:
  void global_refine(int count)
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
      std::vector<double> entityId = generateEntityVisualization(grid_view);
      vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
#if defined(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS) && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
      // boundary id
      std::vector<double> boundaryId = generateBoundaryIdVisualization(grid_view);
      vtkwriter.addCellData(boundaryId, "boundary_id__level_" + Common::to_string(lvl));
#endif
      // write
      vtkwriter.write(filename + "__level_" + Common::to_string(lvl), VTK::appendedraw);
    }
  } // ... visualize_plain(...)

  void visualize_with_boundary(const Common::Configuration& boundary_info_cfg, const std::string filename) const
  {
    if (GridType::dimension > 3) // give us a call if you have any idea!
      DUNE_THROW(NotImplemented, "For grids of dimension > 3!");
    // boundary info
    typedef XT::Grid::BoundaryInfoFactory<typename LevelGridViewType::Intersection> BoundaryInfoFactory;
    auto boundary_info_ptr = BoundaryInfoFactory::create(boundary_info_cfg.get<std::string>("type"), boundary_info_cfg);
    for (auto lvl : Common::value_range(max_level() + 1)) {
      auto grid_view = level_view(lvl);
      // vtk writer
      Dune::VTKWriter<LevelGridViewType> vtkwriter(grid_view);
      // codim 0 entity id
      std::vector<double> entityId = generateEntityVisualization(grid_view);
      vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
#if defined(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS) && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
      // boundary id
      std::vector<double> boundaryId = generateBoundaryIdVisualization(grid_view);
      vtkwriter.addCellData(boundaryId, "boundary_id__level_" + Common::to_string(lvl));
#endif
      // dirichlet values
      std::vector<double> dirichlet = generateBoundaryVisualization(grid_view, *boundary_info_ptr, "dirichlet");
      vtkwriter.addCellData(dirichlet, "isDirichletBoundary__level_" + Common::to_string(lvl));
      // neumann values
      std::vector<double> neumann = generateBoundaryVisualization(grid_view, *boundary_info_ptr, "neumann");
      vtkwriter.addCellData(neumann, "isNeumannBoundary__level_" + Common::to_string(lvl));
      // write
      vtkwriter.write(filename + "__level_" + Common::to_string(lvl), VTK::appendedraw);
    }
  } // ... visualize_with_boundary(...)

#if defined(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS) && DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  std::vector<double> generateBoundaryIdVisualization(const LevelGridViewType& gridView) const
  {
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    const auto it_end = gridView.template end<0>();
    for (auto it = gridView.template begin<0>(); it != it_end; ++it) {
      const auto& entity = *it;
      const auto& index = gridView.indexSet().index(entity);
      data[index] = 0.0;
      size_t numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      const auto intersectionItEnd = gridView.iend(entity);
      for (auto intersectionIt = gridView.ibegin(entity); intersectionIt != intersectionItEnd; ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          data[index] += double(intersectionIt->boundaryId());
        }
      }
      if (isOnBoundary) {
        data[index] /= double(numberOfBoundarySegments);
      }
    } // walk the grid
    return data;
  }
#else
  std::vector<double> generateBoundaryIdVisualization(const LevelGridViewType& /*gridView*/) const
  {
    DUNE_THROW(InvalidStateException, "Experimental Grid Extensions missing");
  }
#endif

  template <class BoundaryInfoType>
  std::vector<double> generateBoundaryVisualization(const LevelGridViewType& gridView,
                                                    const BoundaryInfoType& boundaryInfo,
                                                    const std::string type) const
  {
    constexpr DirichletBoundary dirichlet_type{};
    constexpr NeumannBoundary neumann_type{};
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    for (auto&& entity : elements(gridView)) {
      const auto& index = gridView.indexSet().index(entity);
      data[index] = 0.0;
      for (auto intersectionIt = gridView.ibegin(entity); intersectionIt != gridView.iend(entity); ++intersectionIt) {
        if (type == "dirichlet") {
          if (boundaryInfo.type(*intersectionIt) == dirichlet_type)
            data[index] = 1.0;
        } else if (type == "neumann") {
          if (boundaryInfo.type(*intersectionIt) == neumann_type)
            data[index] = 1.0;
        } else
          DUNE_THROW(Common::Exceptions::internal_error, "Unknown type '" << type << "'!");
      }
    } // walk the grid
    return data;
  } // std::vector< double > generateBoundaryVisualization(...) const

  std::vector<double> generateEntityVisualization(const LevelGridViewType& gridView) const
  {
    std::vector<double> data(gridView.indexSet().size(0));
    for (auto&& entity : elements(gridView)) {
      const auto& index = gridView.indexSet().index(entity);
      data[index] = double(index);
    } // walk the grid
    return data;
  } // ... generateEntityVisualization(...)

  std::shared_ptr<GridType> grid_ptr_;
}; // class GridProvider


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_PROVIDER_HH
