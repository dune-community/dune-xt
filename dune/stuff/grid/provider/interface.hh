
#ifndef DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
#define DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH

// dune-common
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

// dune-grid
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune {

namespace Stuff {

namespace Grid {

namespace Provider {

/**
 *  \brief      Interface for all grid providers.
 *
 *  \attention  The method id() has to be implemented by the derived class, althoug it is not marked virtual (since it
 *              has to be static)
 **/
#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Interface
{
public:
  typedef GridImp GridType;

  typedef Dune::FieldVector<typename GridType::ctype, GridType::dimension> CoordinateType;

  static const std::string id()
  {
    return "stuff.grid.provider.interface";
  }

  virtual GridType& grid() = 0;

  virtual const GridType& grid() const = 0;

  virtual Dune::shared_ptr<GridType> gridPtr() = 0;

  virtual const Dune::shared_ptr<const GridType> gridPtr() const = 0;

private:
  typedef typename GridType::LeafGridView GridViewType;

public:
  void visualize(const std::string filename = id + ".grid") const
  {
    // vtk writer
    GridViewType gridView = grid().leafView();
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
          data[index] += double(intersectionIt->boundaryId());
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

#endif // DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
