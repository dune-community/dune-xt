// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_GRID_ENTITY_VISUALIZATION_HH
#define DUNE_STUFF_GRID_ENTITY_VISUALIZATION_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/geometry/genericgeometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune {
namespace Stuff {
namespace Grid {

struct ElementVisualization
{
  //! Parameter for mapper class
  template <int dim>
  struct P0Layout
  {
    bool contains(const Dune::GeometryType& gt)
    {
      if (gt.dim() == dim)
        return true;
      return false;
    }
  };

  // demonstrate attaching data to elements
  template <class Grid, class F>
  static void elementdata(const Grid& grid, const F& f)
  {
    // get grid view on leaf part
    const auto gridView = grid.leafView();

    // make a mapper for codim 0 entities in the leaf grid
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid, P0Layout> mapper(grid);

    std::vector<double> values(mapper.size());
    for (const auto& entity : DSC::viewRange(gridView)) {
      values[mapper.map(entity)] = f(entity);
    }

    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(gridView);
    vtkwriter.addCellData(values, "data");
    const std::string piecefilesFolderName = "piecefiles";
    const std::string piecefilesPath = f.dir() + "/" + piecefilesFolderName + "/";
    DSC::testCreateDirectory(piecefilesPath);
    vtkwriter.pwrite(f.filename(), f.dir(), piecefilesFolderName, Dune::VTK::appendedraw);
  }


  class FunctorBase
  {
  public:
    FunctorBase(const std::string fname, const std::string dname)
      : filename_(fname)
      , dir_(dname)
    {
    }
    const std::string filename() const
    {
      return filename_;
    }
    const std::string dir() const
    {
      return dir_;
    }

  protected:
    const std::string filename_;
    const std::string dir_;
  };

  class VolumeFunctor : public FunctorBase
  {
  public:
    VolumeFunctor(const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
    {
    }

    template <class Entity>
    double operator()(const Entity& ent) const
    {
      return ent.geometry().volume();
    }
  };

  class ProcessIdFunctor : public FunctorBase
  {
  public:
    ProcessIdFunctor(const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
    {
    }

    template <class Entity>
    double operator()(const Entity& /*ent*/) const
    {
      return Dune::MPIHelper::getCollectiveCommunication().rank();
    }
  };

  template <class GridType>
  class BoundaryFunctor : public FunctorBase
  {
    const GridType& grid_;

  public:
    BoundaryFunctor(const GridType& grid, const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
      , grid_(grid)
    {
    }

    template <class Entity>
    double operator()(const Entity& entity) const
    {
      double ret(0.0);
      int numberOfBoundarySegments(0);
      bool isOnBoundary   = false;
      const auto leafview = grid_.leafView();
      for (const auto& intersection : DSC::intersectionRange(leafview, entity)) {
        if (!intersection.neighbor() && intersection.boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          ret += double(intersection.boundaryId());
        }
      }
      if (isOnBoundary) {
        ret /= double(numberOfBoundarySegments);
      }
      return ret;
    }
  };

  class AreaMarker : public FunctorBase
  {

  public:
    AreaMarker(const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
    {
    }

    template <class Entity>
    double operator()(const Entity& entity) const
    {
      typedef typename Entity::Geometry EntityGeometryType;

      typedef Dune::FieldVector<typename EntityGeometryType::ctype, EntityGeometryType::coorddimension> DomainType;

      const EntityGeometryType& geometry = entity.geometry();

      DomainType baryCenter(0.0);

      for (int corner = 0; corner < geometry.corners(); ++corner) {
        baryCenter += geometry.corner(corner);
      }
      baryCenter /= geometry.corners();

      double ret(0.0);

      if (!((baryCenter[0] < 0.0) || (baryCenter[0] > 1.0))) { // only in unit square
        if (!((baryCenter[1] < 0.0) || (baryCenter[1] > 1.0))) {
          ret = 1.0;
        }
      }
      return ret;
    }
  };

  class GeometryFunctor : public FunctorBase
  {
  public:
    GeometryFunctor(const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
    {
    }

    template <class Entity>
    double operator()(const Entity& ent) const
    {
      const typename Entity::Geometry& geo = ent.geometry();
      double vol = geo.volume();
      if (vol < 0) {
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::setw(8);
        // std::cout.showpoint();
        for (int i = 0; i < geo.corners(); ++i) {
          std::cout << geo.corner(i) << "\t\t";
        }
        std::cout << std::endl;
      }
      return vol;
    }
  };

  //! supply functor
  template <class Grid>
  static void all(const Grid& grid, const std::string outputDir = "visualisation")
  {
    // make function objects
    BoundaryFunctor<Grid> boundaryFunctor(grid, "boundaryFunctor", outputDir);
    AreaMarker areaMarker("areaMarker", outputDir);
    GeometryFunctor geometryFunctor("geometryFunctor", outputDir);
    ProcessIdFunctor processIdFunctor("ProcessIdFunctor", outputDir);
    VolumeFunctor volumeFunctor("volumeFunctor", outputDir);

    // call the visualization functions
    elementdata(grid, boundaryFunctor);
    elementdata(grid, areaMarker);
    elementdata(grid, geometryFunctor);
    elementdata(grid, processIdFunctor);
    elementdata(grid, volumeFunctor);
  }
};


} // namespace Stuff
} // namespace Grid
} // namespace Dune

#endif // DUNE_STUFF_GRID_ENTITY_VISUALIZATION_HH
