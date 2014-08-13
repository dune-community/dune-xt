// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
#define DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH

#include <memory>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_GRID
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

#if HAVE_DUNE_GRID


template <class GridImp>
class ConstProviderInterface
{
public:
  typedef GridImp GridType;
  static const unsigned int dimDomain = GridImp::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  template <ChooseLayer layer, ChoosePartView type>
  struct Layer
  {
    typedef typename Grid::Layer<GridType, layer, type>::Type Type;
  };

  template <ChoosePartView type>
  struct Level
  {
    typedef typename LevelPartView<GridType, type>::Type Type;
  };

  template <ChoosePartView type>
  struct Leaf
  {
    typedef typename LeafPartView<GridType, type>::Type Type;
  };

  template <ChooseLayer type>
  struct View
  {
    typedef typename LayerView<GridType, type>::Type Type;
  };

  template <ChooseLayer type>
  struct Part
  {
    typedef typename LayerView<GridType, type>::Type Type;
  };

  typedef typename Layer<ChooseLayer::level, ChoosePartView::view>::Type LevelGridViewType;
#if HAVE_DUNE_FEM
  typedef typename Layer<ChooseLayer::level, ChoosePartView::part>::Type LevelGridPartType;
#endif

  typedef typename Layer<ChooseLayer::leaf, ChoosePartView::view>::Type LeafGridViewType;
#if HAVE_DUNE_FEM
  typedef typename Layer<ChooseLayer::leaf, ChoosePartView::part>::Type LeafGridPartType;
#endif

  static const std::string static_id()
  {
    return "stuff.grid.provider";
  }

  virtual ~ConstProviderInterface()
  {
  }

  virtual const GridType& grid() const = 0;

  template <ChooseLayer layer_type, ChoosePartView part_view_type>
  std::shared_ptr<const typename Layer<layer_type, part_view_type>::Type> layer(const int level_in = 0) const
  {
    GridType& non_const_grid = const_cast<GridType&>(grid());
    return Grid::Layer<GridType, layer_type, part_view_type>::create(non_const_grid, level_in);
  }

  template <ChoosePartView type>
  std::shared_ptr<const typename Level<type>::Type> level(const int level_in) const
  {
    GridType& non_const_grid = const_cast<GridType&>(grid());
    return LevelPartView<GridType, type>::create(non_const_grid, level_in);
  }

  std::shared_ptr<const LevelGridViewType> level_view(const int level_in) const
  {
    return this->template level<ChoosePartView::view>(level_in);
  }

#if HAVE_DUNE_FEM
  std::shared_ptr<const LevelGridPartType> level_part(const int level_in) const
  {
    return this->template level<ChoosePartView::part>(level_in);
  }
#endif // HAVE_DUNE_FEM

  template <ChoosePartView type>
  std::shared_ptr<const typename Leaf<type>::Type> leaf() const
  {
    GridType& non_const_grid = const_cast<GridType&>(grid());
    return LeafPartView<GridType, type>::create(non_const_grid);
  }

  std::shared_ptr<const LeafGridViewType> leaf_view() const
  {
    return this->template leaf<ChoosePartView::view>();
  }

#if HAVE_DUNE_FEM
  std::shared_ptr<const LeafGridPartType> leaf_part() const
  {
    return this->template leaf<ChoosePartView::part>();
  }
#endif // HAVE_DUNE_FEM

  virtual void visualize(const std::string filename = static_id()) const
  {
    // vtk writer
    auto grid_view = leaf_view();
    Dune::VTKWriter<LeafGridViewType> vtkwriter(*grid_view);
    // codim 0 entity id
    std::vector<double> entityId = generateEntityVisualization(*grid_view);
    vtkwriter.addCellData(entityId, "entity_id");
    // boundary id
    std::vector<double> boundaryId = generateBoundaryIdVisualization(*grid_view);
    vtkwriter.addCellData(boundaryId, "boundary_id");
    // write
    vtkwriter.write(filename, VTK::appendedraw);
  } // ... visualize(...)

  virtual void visualize(const Common::Configuration& boundary_info_cfg, const std::string filename = static_id()) const
  {
    // boundary info
    typedef Stuff::Grid::BoundaryInfoProvider<typename LeafGridViewType::Intersection> BoundaryInfoProvider;
    auto boundary_info_ptr =
        BoundaryInfoProvider::create(boundary_info_cfg.get<std::string>("type"), boundary_info_cfg);
    // vtk writer
    auto grid_view = leaf_view();
    Dune::VTKWriter<LeafGridViewType> vtkwriter(*grid_view);
    // codim 0 entity id
    std::vector<double> entityId = generateEntityVisualization(*grid_view);
    vtkwriter.addCellData(entityId, "entityId");
    // boundary id
    std::vector<double> boundaryId = generateBoundaryIdVisualization(*grid_view);
    vtkwriter.addCellData(boundaryId, "boundaryId");
    // dirichlet values
    std::vector<double> dirichlet = generateBoundaryVisualization(*grid_view, *boundary_info_ptr, "dirichlet");
    vtkwriter.addCellData(dirichlet, "isDirichletBoundary");
    // neumann values
    std::vector<double> neumann = generateBoundaryVisualization(*grid_view, *boundary_info_ptr, "neumann");
    vtkwriter.addCellData(neumann, "isNeumannBoundary");
    // write
    vtkwriter.write(filename, VTK::appendedraw);
  } // ... visualize(...)

private:
  std::vector<double> generateBoundaryIdVisualization(const LeafGridViewType& gridView) const
  {
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    for (typename LeafGridViewType::template Codim<0>::Iterator it = gridView.template begin<0>();
         it != gridView.template end<0>();
         ++it) {
      const auto& entity           = *it;
      const auto& index            = gridView.indexSet().index(entity);
      data[index]                  = 0.0;
      int numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      for (auto intersectionIt = gridView.ibegin(entity); intersectionIt != gridView.iend(entity); ++intersectionIt) {
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
  } // std::vector< double > generateBoundaryIdVisualization(const LeafGridViewType& gridView) const

  template <class BoundaryInfoType>
  std::vector<double> generateBoundaryVisualization(const LeafGridViewType& gridView,
                                                    const BoundaryInfoType& boundaryInfo, const std::string type) const
  {
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    for (const auto& entity : DSC::viewRange(gridView)) {
      const auto& index = gridView.indexSet().index(entity);
      data[index] = 0.0;
      for (auto intersectionIt = gridView.ibegin(entity); intersectionIt != gridView.iend(entity); ++intersectionIt) {
        if (type == "dirichlet") {
          if (boundaryInfo.dirichlet(*intersectionIt))
            data[index] = 1.0;
        } else if (type == "neumann") {
          if (boundaryInfo.neumann(*intersectionIt))
            data[index] = 1.0;
        } else
          DUNE_THROW(Exceptions::internal_error, "Unknown type '" << type << "'!");
      }
    } // walk the grid
    return data;
  } // std::vector< double > generateBoundaryVisualization(...) const

  std::vector<double> generateEntityVisualization(const LeafGridViewType& gridView) const
  {
    std::vector<double> data(gridView.indexSet().size(0));
    for (const auto& entity : DSC::viewRange(gridView)) {
      const auto& index = gridView.indexSet().index(entity);
      data[index]       = double(index);
    } // walk the grid
    return data;
  } // ... generateEntityVisualization(...)
}; // class ConstProviderInterface


template <class GridImp>
class ProviderInterface : public ConstProviderInterface<GridImp>
{
  typedef ConstProviderInterface<GridImp> BaseType;

public:
  using typename BaseType::GridType;
  using typename BaseType::LevelGridViewType;
  using typename BaseType::LevelGridPartType;
  using typename BaseType::LeafGridViewType;
  using typename BaseType::LeafGridPartType;

  static const std::string static_id()
  {
    return BaseType::static_id();
  }

  virtual ~ProviderInterface()
  {
  }

  using BaseType::grid;

  virtual GridType& grid() = 0;

  using BaseType::layer;

  template <ChooseLayer layer_type, ChoosePartView part_view_type>
  std::shared_ptr<typename BaseType::template Layer<layer_type, part_view_type>::Type> layer(const int level_in = 0)
  {
    return Grid::Layer<GridType, layer_type, part_view_type>::create(grid(), level_in);
  }

  using BaseType::level;

  template <ChoosePartView type>
  std::shared_ptr<typename BaseType::template Level<type>::Type> level(const int level_in)
  {
    return LevelPartView<GridType, type>::create(grid(), level_in);
  }

  using BaseType::level_view;

  std::shared_ptr<LevelGridViewType> level_view(const int level_in)
  {
    return this->template level<ChoosePartView::view>(level_in);
  }

#if HAVE_DUNE_FEM
  using BaseType::level_part;

  std::shared_ptr<LevelGridPartType> level_part(const int level_in)
  {
    return this->template level<ChoosePartView::part>(level_in);
  }
#endif // HAVE_DUNE_FEM

  using BaseType::leaf;

  template <ChoosePartView type>
  std::shared_ptr<typename BaseType::template Leaf<type>::Type> leaf()
  {
    return LeafPartView<GridType, type>::create(grid());
  }

  using BaseType::leaf_view;

  std::shared_ptr<LeafGridViewType> leaf_view()
  {
    return this->template leaf<ChoosePartView::view>();
  }

#if HAVE_DUNE_FEM
  using BaseType::leaf_part;

  std::shared_ptr<LeafGridPartType> leaf_part()
  {
    return this->template leaf<ChoosePartView::part>();
  }

#endif // HAVE_DUNE_FEM
}; // class ProviderInterface


#else // HAVE_DUNE_GRID


template <class GridImp>
class ConstProviderInterface
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};


template <class GridImp>
class ProviderInterface
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-grid!");
};


#endif // HAVE_DUNE_GRID

} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
