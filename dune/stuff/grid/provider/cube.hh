#ifndef DUNE_STUFF_GRID_PROVIDER_CUBE_HH
#define DUNE_STUFF_GRID_PROVIDER_CUBE_HH

// system
#include <sstream>

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-grid
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune {

namespace Stuff {

namespace Grid {

namespace Provider {

/**
  \brief  Creates a grid of a cube in various dimensions.

          Default implementation using the Dune::StructuredGridFactory to create a grid of a cube in 1, 2 or 3
          dimensions. Tested with
          <ul><li> \c YASPGRID, \c variant 1, dim = 1, 2, 3,
          <li> \c SGRID, \c variant 1, dim = 1, 2, 3,
          <li> \c ALUGRID_SIMPLEX, \c variant 2, dim = 2, 3,
          <li> \c ALUGRID_CONFORM, \c variant 2, dim = 2, 2 and
          <li> \c ALUGRID_CUBE, \c variant 1, dim = 2, 3.</ul>
  \tparam GridImp
          Type of the underlying grid.
  \tparam variant
          Type of the codim 0 elements:
          <ul><li>\c 1: cubes
          <li>2: simplices</ul>
  **/
template <typename GridImp, int variant>
class GenericCube
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  //! Dimension of the provided grid.
  static const unsigned int dim = GridType::dimension;

  //! Type of the grids coordinates.
  typedef Dune::FieldVector<typename GridType::ctype, dim> CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.cube
  static const std::string id;

  /**
    \brief      Creates a cube.
    \param[in]  paramTree
                A Dune::ParameterTree containing
                <ul><li> the following keys directly or
                <li> a subtree named Cube::id, containing the following keys.</ul>
                The actual keys are:
                <ul><li> \c lowerLeft: \a double that is used as a lower left corner in each dimension.
                <li> \c upperRight: \a double that is used as a upper right corner in each dimension.</ul>
    **/
  GenericCube(const Dune::ParameterTree& paramTree)
    : lowerLeft_(0.0)
    , upperRight_(0.0)
  {
    // get outer bounds
    const double lowerLeft  = paramTree.get("lowerLeft", 0.0);
    const double upperRight = paramTree.get("upperRight", 1.0);
    assert(lowerLeft < upperRight);
    lowerLeft_  = lowerLeft;
    upperRight_ = upperRight;
    // get number of elements per dim
    if (paramTree.hasKey("level"))
      std::fill(numElements_.begin(), numElements_.end(), std::pow(2, paramTree.get("level", 1)));
    else {
      for (unsigned int d = 0; d < dim; ++d) {
        std::stringstream s;
        s << "numElements." << d;
        numElements_[d] = paramTree.get(s.str(), 1);
      }
    }
    buildGrid();
  } // Cube(const Dune::ParameterTree& paramTree)

  /**
    \brief      Creates a cube.
    \param[in]  lowerLeft
                A vector that is used as a lower left corner.
    \param[in]  upperRight
                A vector that is used as a upper right corner.
    **/
  GenericCube(const CoordinateType& lowerLeft, const CoordinateType& upperRight, const int level = 1)
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
  {
    std::fill(numElements_.begin(), numElements_.end(), std::pow(2, level));
    buildGrid();
  }

  /**
    \brief      Creates a cube.
    \param[in]  lowerLeft
                A double that is used as a lower left corner in each dimension.
    \param[in]  upperRight
                A double that is used as a upper right corner in each dimension.
    **/
  GenericCube(const double lowerLeft, const double upperRight, const int level = 1)
    : lowerLeft_(lowerLeft)
    , upperRight_(upperRight)
  {
    std::fill(numElements_.begin(), numElements_.end(), std::pow(2, level));
    buildGrid();
  }

  /**
    \brief  Provides access to the created grid.
    \return Reference to the grid.
    **/
  GridType& grid()
  {
    return *grid_;
  }

  const GridType& grid() const
  {
    return *grid_;
  }

  const CoordinateType& lowerLeft() const
  {
    return lowerLeft_;
  }

  const CoordinateType& upperRight() const
  {
    return upperRight_;
  }

private:
  template <int dim>
  struct P0Layout
  {
    bool contains(Dune::GeometryType& geometry)
    {
      if (geometry.dim() == dim)
        return true;
      return false;
    }
  }; // layout class for codim 0 mapper

  /**
    \brief      Visualizes the grid using Dune::VTKWriter.
    \param[in]  paramTree
                A Dune::ParameterTree containing
                <ul><li> the following keys directly or
                <li> a subtree named Cube::id, containing the following keys, or
                <li> a subtree named Cube::id + \c .visualize, containing the following keys.</ul>
                The actual keys are:
                <ul><li> \c grid: if specified, filename of the vtk file in which the grid which can be obtained via
                  grid() is visualized (\a if \a not \a specified: \a no \a visualization).
                <li> \c mdGrid: if specified, filename of the vtk file in which the multidomain grid which can be
                  obtained via mdGrid() is visualized (\a if \a not \a specified: \a no \a visualization).</ul>
    **/
public:
  void visualize(Dune::ParameterTree& paramTree) const
  {
    const std::string localId = "visualize";
    // get correct parameters
    Dune::ParameterTree* paramsP = &paramTree;
    if (paramsP->hasSub(id))
      paramsP = &(paramsP->sub(id));
    if (paramsP->hasSub(localId))
      paramsP                   = &(paramsP->sub(localId));
    Dune::ParameterTree& params = *paramsP;
    // check for grid visualization
    if (params.hasKey("grid")) {
      const std::string filenameGrid = params.get("grid", id + ".grid");
      // grid view
      typedef GridImp GridType;
      GridType& grid = this->grid();
      typedef typename GridType::LeafGridView GridView;
      GridView gridView = grid.leafView();
      // mapper
      Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType, P0Layout> mapper(grid);
      std::vector<double> data(mapper.size());
      // walk the grid
      typedef typename GridView::template Codim<0>::Iterator ElementIterator;
      typedef typename GridView::template Codim<0>::Entity ElementType;
      typedef typename ElementType::LeafIntersectionIterator FacetIteratorType;
      for (ElementIterator it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it) {
        ElementType& element = *it;
        data[mapper.map(element)] = 0.0;
        int numberOfBoundarySegments = 0;
        bool isOnBoundary = false;
        for (FacetIteratorType facet = element.ileafbegin(); facet != element.ileafend(); ++facet) {
          if (!facet->neighbor() && facet->boundary()) {
            isOnBoundary = true;
            numberOfBoundarySegments += 1;
            data[mapper.map(element)] += double(facet->boundaryId());
          }
        }
        if (isOnBoundary) {
          data[mapper.map(element)] /= double(numberOfBoundarySegments);
        }
      } // walk the grid
      // write to vtk
      Dune::VTKWriter<GridView> vtkwriter(gridView);
      vtkwriter.addCellData(data, "boundaryId");
      vtkwriter.write(filenameGrid, Dune::VTK::ascii);
    } // check for grid visualization
  } // void visualize(Dune::ParameterTree& paramTree)

private:
  void buildGrid()
  {
    switch (variant) {
      case 1:
        grid_ = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerLeft_, upperRight_, numElements_);
        break;
      case 2:
        grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(lowerLeft_, upperRight_, numElements_);
        break;
      default:
        DUNE_THROW(Dune::NotImplemented, "Variant " << variant << " of cube not implemented.");
    }
    return;
  } // void buildGrid(const CoordinateType& lowerLeft, const CoordinateType& upperRight)

  CoordinateType lowerLeft_;
  CoordinateType upperRight_;
  Dune::array<unsigned int, dim> numElements_;
  Dune::shared_ptr<GridType> grid_;
}; // class GenericCube

template <typename GridImp, int variant>
const std::string GenericCube<GridImp, variant>::id = "stuff.grid.provider.cube";

// default implementation of a cube for any grid
// tested for
// dim = 2
//  ALUGRID_SIMPLEX, variant 2
//  ALUGRID_CONFORM, variant 2
// dim = 3
//  ALUGRID_SIMPLEX, variant 2
template <typename GridType>
class Cube : public GenericCube<GridType, 2>
{
private:
  typedef GenericCube<GridType, 2> BaseType;

public:
  typedef typename BaseType::CoordinateType CoordinateType;

  Cube(const Dune::ParameterTree& paramTree)
    : BaseType(paramTree)
  {
  }

  Cube(const CoordinateType& lowerLeft, const CoordinateType& upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }

  Cube(const double lowerLeft, const double upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }
}; // class Cube

// specialization of Cube for YaspGrid
// tested for dim = 1, 2, 3
template <int dim>
class Cube<Dune::YaspGrid<dim>> : public GenericCube<Dune::YaspGrid<dim>, 1>
{
private:
  typedef GenericCube<Dune::YaspGrid<dim>, 1> BaseType;

public:
  typedef typename BaseType::CoordinateType CoordinateType;

  Cube(const Dune::ParameterTree& paramTree)
    : BaseType(paramTree)
  {
  }

  Cube(const CoordinateType& lowerLeft, const CoordinateType& upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }

  Cube(const double lowerLeft, const double upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }
}; // class Cube< Dune::YaspGrid< dim > >

// specialization of Cube for SGrid
// tested for dim = 1, 2, 3
template <int dim>
class Cube<Dune::SGrid<dim, dim>> : public GenericCube<Dune::SGrid<dim, dim>, 1>
{
private:
  typedef GenericCube<Dune::SGrid<dim, dim>, 1> BaseType;

public:
  typedef typename BaseType::CoordinateType CoordinateType;

  Cube(const Dune::ParameterTree& paramTree)
    : BaseType(paramTree)
  {
  }

  Cube(const CoordinateType& lowerLeft, const CoordinateType& upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }

  Cube(const double lowerLeft, const double upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }
}; // class Cube< Dune::SGrid< dim, dim > >

// specialization of Cube for ALUCubeGrid
// tested for dim = 2, 3
template <int dim>
class Cube<Dune::ALUCubeGrid<dim, dim>> : public GenericCube<Dune::ALUCubeGrid<dim, dim>, 1>
{
private:
  typedef GenericCube<Dune::ALUCubeGrid<dim, dim>, 1> BaseType;

public:
  typedef typename BaseType::CoordinateType CoordinateType;

  Cube(const Dune::ParameterTree& paramTree)
    : BaseType(paramTree)
  {
  }

  Cube(const CoordinateType& lowerLeft, const CoordinateType& upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }

  Cube(const double lowerLeft, const double upperRight, const int level = 1)
    : BaseType(lowerLeft, upperRight, level)
  {
  }
}; // class Cube< Dune::ALUCubeGrid< dim, dim > >

template <typename GridType>
class UnitCube : public Cube<GridType>
{
private:
  typedef Cube<GridType> BaseType;

public:
  UnitCube(const Dune::ParameterTree& paramTree)
    : BaseType(0.0, 1.0, paramTree.get("level", 1))
  {
  }

  UnitCube(const int level = 1)
    : BaseType(0.0, 1.0, level)
  {
  }
}; // class UnitCube

} // namespace Provider

} // namespace Grid

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_CUBE_HH
