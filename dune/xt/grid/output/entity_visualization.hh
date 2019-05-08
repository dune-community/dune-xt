// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2013 - 2016, 2018)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_XT_GRID_OUTPUT_ENTITY_VISUALIZATION_HH
#define DUNE_XT_GRID_OUTPUT_ENTITY_VISUALIZATION_HH

#include <boost/io/ios_state.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/ranges.hh>

namespace Dune {
namespace XT {
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
    const auto gridView = grid.leafGridView();

    // make a mapper for codim 0 entities in the leaf grid
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid, P0Layout> mapper(grid);

    std::vector<double> values(mapper.size());
    for (auto&& entity : elements(gridView)) {
      values[mapper.index(entity)] = f(entity);
    }

    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(gridView);
    vtkwriter.addCellData(values, "data");
    const std::string piecefilesFolderName = "piecefiles";
    const std::string piecefilesPath = f.dir() + "/" + piecefilesFolderName + "/";
    Common::test_create_directory(piecefilesPath);
    vtkwriter.pwrite(f.filename(), f.dir(), piecefilesFolderName, Dune::VTK::appendedraw);
  }

  class FunctorBase
  {
  public:
    FunctorBase(const std::string fname, const std::string dname)
      : filename_(fname)
      , dir_(dname)
    {}
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
    {}

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
    {}

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
    {}

    template <class Entity>
    double operator()(const Entity& entity) const
    {
      double ret(0.0);
      int numberOfBoundarySegments(0);
      bool isOnBoundary = false;
      const auto leafview = grid_.leafGridView();
      const auto intersection_it_end = leafview.iend(entity);
      for (auto intersection_it = leafview.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
        const auto& intersection = *intersection_it;
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
    {}

    template <class Entity>
    double operator()(const Entity& entity) const
    {
      typedef typename Entity::Geometry EntityGeometryType;

      typedef Dune::FieldVector<typename EntityGeometryType::ctype, EntityGeometryType::coorddimension> DomainType;

      const EntityGeometryType& geometry = entity.geometry();

      DomainType baryCenter(0.0);

      for (auto corner : Common::value_range(geometry.corners())) {
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
    {}

    template <class Entity>
    double operator()(const Entity& ent) const
    {
      const typename Entity::Geometry& geo = ent.geometry();
      double vol = geo.volume();
      if (vol < 0) {
        boost::io::ios_all_saver guard(DXTC_LOG_ERROR);
        DXTC_LOG_ERROR << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::setw(8);
        // std::cout.showpoint();
        for (auto i : Common::value_range(geo.corners())) {
          DXTC_LOG_ERROR << geo.corner(i) << "\t\t";
        }
        DXTC_LOG_ERROR << std::endl;
      }
      return vol;
    }
  };

  class PartitionTypeFunctor : public FunctorBase
  {
  public:
    PartitionTypeFunctor(const std::string fname, const std::string dname)
      : FunctorBase(fname, dname)
    {}

    template <class Entity>
    double operator()(const Entity& ent) const
    {
      const typename Entity::Geometry& geo = ent.geometry();
      const int type{static_cast<int>(ent.partitionType())};
      DXTC_LOG_ERROR << "TYPE " << type << std::endl;
      return static_cast<double>(type);
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
    PartitionTypeFunctor partitionTypeFunctor("partitionTypeFunctor", outputDir);

    // call the visualization functions
    elementdata(grid, boundaryFunctor);
    elementdata(grid, areaMarker);
    elementdata(grid, geometryFunctor);
    elementdata(grid, processIdFunctor);
    elementdata(grid, volumeFunctor);
    elementdata(grid, partitionTypeFunctor);
  }
};

template <class GridType>
void visualize_index_per_level(const GridType& grid_, std::string filename)
{
  if (GridType::dimension > 3)
    DUNE_THROW(NotImplemented, "For grids of dimension > 3!");
  for (auto lvl : Common::value_range(grid_.maxLevel() + 1)) {
    const auto grid_view = grid_.levelGridView(lvl);
    std::vector<double> entityId(grid_view.indexSet().size(0));
    for (auto&& entity : elements(grid_view)) {
      const auto& index = grid_view.indexSet().index(entity);
      entityId[index] = double(index);
    }
    Dune::VTKWriter<decltype(grid_view)> vtkwriter(grid_view);
    vtkwriter.addCellData(entityId, "entity_id__level_" + Common::to_string(lvl));
    vtkwriter.write(filename + "__level_" + Common::to_string(lvl), VTK::appendedraw);
  }
} // ... visualize_plain(...)

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_OUTPUT_ENTITY_VISUALIZATION_HH
