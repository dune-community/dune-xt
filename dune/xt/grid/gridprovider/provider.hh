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
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_XT_GRID_PROVIDER_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_PROVIDER_HH

#include <memory>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/capabilities.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {

#if HAVE_DUNE_FEM
template <class T>
using DefaultDDGridImp = DD::SubdomainGrid<T>;
#else
template <class T>
using DefaultDDGridImp = AlwaysFalse<T>;
#endif

template <class GridImp, typename DdGridImp = DefaultDDGridImp<GridImp>>
class GridProvider
{
  static_assert(is_grid<GridImp>::value, "");

public:
  typedef GridProvider<GridImp, DdGridImp> ThisType;

  typedef GridImp GridType;
  typedef DdGridImp DdGridType;
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
  GridProvider(GridType* grd_ptr, DdGridType* dd_grd_ptr = nullptr)
    : grid_ptr_(grd_ptr)
    , dd_grid_ptr_(dd_grd_ptr)
  {
  }

  GridProvider(std::shared_ptr<GridType> grd_ptr, std::shared_ptr<DdGridType> dd_grd_ptr = nullptr)
    : grid_ptr_(grd_ptr)
    , dd_grid_ptr_(dd_grd_ptr)
  {
  }

  GridProvider(const ThisType& other) = default;

  // Manual ctor required for clang 3.8.1-12~bpo8+1 (otherwise: undefined reference).
  GridProvider(ThisType&& source)
    : grid_ptr_(source.grid_ptr_)
    , dd_grid_ptr_(source.dd_grid_ptr_)
  {
  }

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

  const std::shared_ptr<DdGridType>& dd_grid_ptr() const
  {
    if (!dd_grid_ptr_)
      DUNE_THROW(InvalidStateException, "No DD grid provided on construction!");
    return dd_grid_ptr_;
  }

  std::shared_ptr<DdGridType> dd_grid_ptr()
  {
    if (!dd_grid_ptr_)
      DUNE_THROW(InvalidStateException, "No DD grid provided on construction!");
    return dd_grid_ptr_;
  }

  const DdGridType& dd_grid() const
  {
    if (!dd_grid_ptr_)
      DUNE_THROW(InvalidStateException, "No DD grid provided on construction!");
    return *dd_grid_ptr_;
  }

  DdGridType& dd_grid()
  {
    if (!dd_grid_ptr_)
      DUNE_THROW(InvalidStateException, "No DD grid provided on construction!");
    return *dd_grid_ptr_;
  }

  int max_level() const
  {
    return grid_ptr_->maxLevel();
  }

  template <Backends backend>
  typename Layer<GridType, Layers::level, backend>::type DUNE_DEPRECATED_MSG("Use level_view() instead (09.05.2017)!")
      level(const int lvl = 0) const
  {
    return Layer<GridType, Layers::level, backend>::create(*grid_ptr_, lvl);
  }

  template <Backends backend>
  typename Layer<GridType, Layers::leaf, backend>::type DUNE_DEPRECATED_MSG("Use leaf_view() instead (09.05.2017)!")
      leaf() const
  {
    return Layer<GridType, Layers::leaf, backend>::create(*grid_ptr_);
  }

  template <Layers lr, Backends backend>
  typename Layer<GridType, lr, backend, DdGridType>::type layer(const int lvl = 0) const
  {
    return Layer<GridType, lr, backend, DdGridType>::create(*grid_ptr_, lvl, dd_grid_ptr_);
  }

  template <Layers lr>
  typename Layer<GridType, lr, Backends::part, DdGridType>::type layer(const int lvl = 0)
  {
    return Layer<GridType, lr, Backends::part, DdGridType>::create(*grid_ptr_, lvl, dd_grid_ptr_);
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
  template <bool is_dd_subdomain = std::is_same<DdGridType, DD::SubdomainGrid<GridType>>::value, bool anything = true>
  struct visualize_dd_helper;

#if HAVE_DUNE_FEM
  template <bool anything>
  struct visualize_dd_helper<true, anything>
  {
    void
    operator()(const std::shared_ptr<DdGridType>& dd_grid_ptr, const std::string filename, const bool with_coupling)
    {
      if (!dd_grid_ptr)
        DUNE_THROW(InvalidStateException, "No DD grid provided on construction!");
      const auto& dd_grid = *dd_grid_ptr;
      // vtk writer
      typedef typename DdGridType::GlobalGridPartType GlobalGridPartType;
      const auto& globalGridPart = dd_grid.globalGridPart();
      typedef Dune::Fem::GridPart2GridView<GlobalGridPartType> GVT;
      GVT globalGridView(globalGridPart);
      Dune::VTKWriter<GVT> vtkwriter(globalGridView);
      // data
      std::map<std::string, std::vector<double>> data;
      data["subdomain"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
      data["global entity id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
      data["global boundary id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
      data["local boundary id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
      // walk the global grid view
      for (auto it = globalGridView.template begin<0>(); it != globalGridView.template end<0>(); ++it) {
        const auto& entity = *it;
        const auto index = globalGridView.indexSet().index(entity);
        data["subdomain"][index] = dd_grid.subdomainOf(index);
        data["global entity id"][index] = double(index);
        // compute global boundary id
        data["global boundary id"][index] = 0.0;
        int numberOfBoundarySegments = 0;
        bool isOnBoundary = false;
        for (auto intersectionIt = globalGridView.ibegin(entity); intersectionIt != globalGridView.iend(entity);
             ++intersectionIt) {
          if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
            isOnBoundary = true;
            numberOfBoundarySegments += 1;
            data["global boundary id"][index] += double(intersectionIt->boundarySegmentIndex());
          }
        }
        if (isOnBoundary) {
          data["global boundary id"][index] /= double(numberOfBoundarySegments);
        } // compute global boundary id
      } // walk the global grid view
      // walk the subdomains
      for (unsigned int s = 0; s < dd_grid.size(); ++s) {
        // walk the local grid view
        const auto localGridView = dd_grid.localGridPart(s);
        for (auto it = localGridView.template begin<0>(); it != localGridView.template end<0>(); ++it) {
          const auto& entity = *it;
          const unsigned int index = globalGridView.indexSet().index(entity);
          // compute local boundary id
          unsigned int numberOfBoundarySegments = 0u;
          for (auto intersectionIt = localGridView.ibegin(entity); intersectionIt != localGridView.iend(entity);
               ++intersectionIt) {
            if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
              numberOfBoundarySegments += 1u;
              data["local boundary id"][index] += double(intersectionIt->boundarySegmentIndex());
            }
          }
          if (numberOfBoundarySegments > 0)
            data["local boundary id"][index] /= double(numberOfBoundarySegments);
          // visualize coupling
          if (with_coupling) {
            for (auto nn : dd_grid.neighborsOf(s)) {
              const auto coupling_grid_view = dd_grid.couplingGridPart(s, nn);
              const std::string coupling_str = "coupling (" + Common::to_string(s) + ", " + Common::to_string(nn) + ")";
              data[coupling_str] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
              const auto entity_it_end = coupling_grid_view.template end<0>();
              for (auto entity_it = coupling_grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
                const auto& coupling_entity = *entity_it;
                const size_t global_entity_id = globalGridView.indexSet().index(coupling_entity);
                data[coupling_str][global_entity_id] = 1.0;
                for (auto intersection_it = coupling_grid_view.ibegin(coupling_entity);
                     intersection_it != coupling_grid_view.iend(coupling_entity);
                     ++intersection_it) {
                  const auto& intersection = *intersection_it;
                  if (intersection.neighbor() && !intersection.boundary()) {
                    const auto neighbor = intersection.outside();
                    const size_t global_neighbor_id = globalGridView.indexSet().index(neighbor);
                    data[coupling_str][global_neighbor_id] = 0.5;
                  }
                }
              }
            }
          } // if (with_coupling)
        } // walk the local grid view
      } // walk the subdomains
      if (dd_grid.oversampling()) {
        // walk the oversampled grid parts
        for (size_t ss = 0; ss < dd_grid.size(); ++ss) {
          const std::string string_id = "oversampled subdomain " + Common::to_string(ss);
          data[string_id] = std::vector<double>(globalGridView.indexSet().size(0), -1.0);
          typedef typename DdGridType::LocalGridPartType LocalGridPartType;
          const LocalGridPartType oversampledGridPart = dd_grid.localGridPart(ss, true);
          for (auto it = oversampledGridPart.template begin<0>(); it != oversampledGridPart.template end<0>(); ++it) {
            const auto& entity = *it;
            const auto index = globalGridView.indexSet().index(entity);
            data[string_id][index] = 0.0;
            // compute local boundary id
            unsigned int numberOfBoundarySegments = 0u;
            for (auto intersectionIt = oversampledGridPart.ibegin(entity);
                 intersectionIt != oversampledGridPart.iend(entity);
                 ++intersectionIt) {
              if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
                numberOfBoundarySegments += 1u;
                data[string_id][index] += double(intersectionIt->boundarySegmentIndex());
              }
            }
            if (numberOfBoundarySegments > 0) {
              data[string_id][index] /= double(numberOfBoundarySegments);
            } // compute global boundary id
          }
        }
      } // if (dd_grid.oversampling())

      // write
      for (const auto& data_pair : data)
        vtkwriter.addCellData(data_pair.second, data_pair.first);
      vtkwriter.write(filename, Dune::VTK::appendedraw);
    } // ... operator()(...)
  }; // struct visualize_dd_helper<true, ...>
#endif // HAVE_DUNE_FEM

  template <bool anything>
  struct visualize_dd_helper<false, anything>
  {
    void operator()(const std::shared_ptr<DdGridType>& /*dd_grid_ptr*/,
                    const std::string /*filename*/,
                    const bool /*with_coupling*/)
    {
      DUNE_THROW(NotImplemented, "Only available for DD::SubdomainGrid!");
    }
  };

public:
  void visualize_dd(const std::string filename, const bool with_coupling) const
  {
    visualize_dd_helper<>()(dd_grid_ptr_, filename, with_coupling);
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
  void global_refine(const int count = 1)
  {
    global_refine_helper<GridType>()(grid(), count);
  }

private:
  template <class G, bool enable = has_boundary_id<G>::value>
  struct add_boundary_id_visualization
  {
    add_boundary_id_visualization()
    {
    }

    template <class V>
    void operator()(V& vtk_writer, const std::vector<double>& boundary_id, const int lvl) const
    {
      vtk_writer.addCellData(boundary_id, "boundary_id__level_" + Common::to_string(lvl));
    }

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
  }; // struct add_boundary_id_visualization<..., true>

  template <class G>
  struct add_boundary_id_visualization<G, false>
  {
    add_boundary_id_visualization()
    {
    }

    template <class V>
    void operator()(V& /*vtk_writer*/, const std::vector<double>& /*boundary_id*/, const int /*lvl*/) const
    {
    }

    std::vector<double> generateBoundaryIdVisualization(const LevelGridViewType&) const
    {
      return std::vector<double>();
    }
  };

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
      // boundary id
      const add_boundary_id_visualization<GridType> add_boundary_id;
      const std::vector<double> boundary_id = add_boundary_id.generateBoundaryIdVisualization(grid_view);
      add_boundary_id(vtkwriter, boundary_id, lvl);
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
    auto boundary_info_ptr = BoundaryInfoFactory::create(boundary_info_cfg);
    for (auto lvl : Common::value_range(max_level() + 1)) {
      auto grid_view = level_view(lvl);
      // vtk writer
      Dune::VTKWriter<LevelGridViewType> vtkwriter(grid_view);
      // codim 0 entity id
      std::vector<double> entityId = generateEntityVisualization(grid_view);
      vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
      // boundary id
      const add_boundary_id_visualization<GridType> add_boundary_id;
      const std::vector<double> boundary_id = add_boundary_id.generateBoundaryIdVisualization(grid_view);
      add_boundary_id(vtkwriter, boundary_id, lvl);
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
  std::shared_ptr<DdGridType> dd_grid_ptr_;
}; // class GridProvider


} // namespace Grid
} // namespace XT
} // namespace Dune


#include "provider.lib.hh"


#endif // DUNE_XT_GRID_PROVIDER_PROVIDER_HH
