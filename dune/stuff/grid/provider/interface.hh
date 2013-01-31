#ifndef DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
#define DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/parametertree.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Provider {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Interface
{
public:
  typedef GridImp GridType;

  typedef Interface<GridType> ThisType;

  static const unsigned int dim = GridType::dimension;

  typedef typename GridType::ctype ctype;

  typedef Dune::FieldVector<ctype, dim> CoordinateType;

  static const std::string id()
  {
    return "stuff.grid.provider";
  }

  virtual Dune::shared_ptr<GridType> grid() = 0;

  virtual const Dune::shared_ptr<const GridType> grid() const = 0;

private:
  typedef typename GridType::LeafGridView GridViewType;

public:
  virtual void visualize(const std::string filename = id()) const
  {
    // vtk writer
    GridViewType gridView = grid()->leafView();
    Dune::VTKWriter<GridViewType> vtkwriter(gridView);
    // boundary id
    std::vector<double> boundaryId = generateBoundaryIdVisualization(gridView);
    vtkwriter.addCellData(boundaryId, "boundaryId");
    // codim 0 entity id
    std::vector<double> entityId = generateEntityVisualization(gridView);
    vtkwriter.addCellData(entityId, "entityId");
    // write
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const std::string filename = id + ".grid") const

private:
  std::vector<double> generateBoundaryIdVisualization(const GridViewType& gridView) const
  {
    typedef typename GridViewType::IndexSet::IndexType IndexType;
    typedef typename GridViewType::template Codim<0>::Entity EntityType;
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    for (typename GridViewType::template Codim<0>::Iterator it = gridView.template begin<0>();
         it != gridView.template end<0>();
         ++it) {
      const EntityType& entity     = *it;
      const IndexType& index       = gridView.indexSet().index(entity);
      data[index]                  = 0.0;
      int numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      for (typename GridViewType::IntersectionIterator intersectionIt = gridView.ibegin(entity);
           intersectionIt != gridView.iend(entity);
           ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
          data[index] += double(intersectionIt->boundaryId());
#pragma GCC diagnostic pop
        }
      }
      if (isOnBoundary) {
        data[index] /= double(numberOfBoundarySegments);
      }
    } // walk the grid
    return data;
  } // std::vector< double > generateBoundaryIdVisualization(const GridViewType& gridView) const

  std::vector<double> generateEntityVisualization(const GridViewType& gridView) const
  {
    typedef typename GridViewType::IndexSet::IndexType IndexType;
    typedef typename GridViewType::template Codim<0>::Entity EntityType;
    std::vector<double> data(gridView.indexSet().size(0));
    // walk the grid
    for (typename GridViewType::template Codim<0>::Iterator it = gridView.template begin<0>();
         it != gridView.template end<0>();
         ++it) {
      const EntityType& entity = *it;
      const IndexType& index   = gridView.indexSet().index(entity);
      data[index]              = double(index);
    } // walk the grid
    return data;
  } // std::vector< double > generateEntityVisualization(const GridViewType& gridView) const
}; // class Interface

} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
