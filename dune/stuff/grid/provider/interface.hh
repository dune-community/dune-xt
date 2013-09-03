#ifndef DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
#define DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH

#include <config.h>
#if HAVE_DUNE_GRID

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/ranges.hh>

namespace Dune {
namespace Stuff {


template <class GridImp = Dune::SGrid<2, 2>>
class GridProviderInterface
{
public:
  static const unsigned int dim = GridImp::dimension;

  typedef GridImp GridType;
  typedef GridProviderInterface<GridType> ThisType;
  typedef typename GridType::ctype ctype;
  typedef Dune::FieldVector<ctype, dim> CoordinateType;

private:
  typedef typename GridType::LeafGridView LeafGridViewType;

public:
  static const std::string id()
  {
    return "gridprovider";
  }

  virtual ~GridProviderInterface()
  {
  }

  virtual std::shared_ptr<GridType> grid() = 0;
  virtual const std::shared_ptr<const GridType> grid() const = 0;

  virtual void visualize(const std::string filename = id()) const
  {
    // vtk writer
    LeafGridViewType gridView = grid()->leafView();
    Dune::VTKWriter<LeafGridViewType> vtkwriter(gridView);
    // boundary id
    std::vector<double> boundaryId = generateBoundaryIdVisualization(gridView);
    vtkwriter.addCellData(boundaryId, "boundaryId");
    // codim 0 entity id
    std::vector<double> entityId = generateEntityVisualization(gridView);
    vtkwriter.addCellData(entityId, "entityId");
    // write
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const std::string filename = id + ".grid") const

  virtual void visualize(const std::string boundaryInfoType, const Dune::ParameterTree& description,
                         const std::string filename = id()) const
  {
    // vtk writer
    LeafGridViewType gridView = grid()->leafView();
    Dune::VTKWriter<LeafGridViewType> vtkwriter(gridView);
    // boundary id
    std::vector<double> boundaryId = generateBoundaryIdVisualization(gridView);
    vtkwriter.addCellData(boundaryId, "boundaryId");
    const GridboundaryInterface<typename LeafGridViewType::Intersection>* boundaryInfo =
        Gridboundaries<typename LeafGridViewType::Intersection>::create(boundaryInfoType, description);
    // dirichlet values
    std::vector<double> dirichlet = generateBoundaryVisualization(gridView, *boundaryInfo, "dirichlet");
    vtkwriter.addCellData(dirichlet, "isDirichletBoundary");
    // neumann values
    std::vector<double> neumann = generateBoundaryVisualization(gridView, *boundaryInfo, "neumann");
    delete boundaryInfo;
    vtkwriter.addCellData(neumann, "isNeumannBoundary");
    // codim 0 entity id
    std::vector<double> entityId = generateEntityVisualization(gridView);
    vtkwriter.addCellData(entityId, "entityId");
    // write
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const std::string filename = id + ".grid") const

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

  std::vector<double> generateBoundaryVisualization(
      const LeafGridViewType& gridView,
      const Dune::Stuff::GridboundaryInterface<typename LeafGridViewType::Intersection>& boundaryInfo,
      const std::string type) const
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
          DUNE_THROW(Dune::InvalidStateException, "BOOM");
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
  } // std::vector< double > generateEntityVisualization(const LeafGridViewType& gridView) const
}; // class GridProviderInterface


} // namespace Provider
} // namespace Grid

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
