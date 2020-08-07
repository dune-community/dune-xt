// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017, 2019)
//   René Fritze     (2013 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2016, 2019 - 2020)

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
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/boundaryinfo/types.hh>
#include <dune/xt/grid/capabilities.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune::XT::Grid {

struct ElementVisualization
{
  // demonstrate attaching data to elements
  template <class View, class F>
  static void elementdata(const View& view, const F& f)
  {
    // make a mapper for codim 0 entities in the leaf grid
    using Grid = extract_grid_t<View>;
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid> mapper(view.grid(), mcmgElementLayout());

    std::vector<double> values(mapper.size());
    for (auto&& entity : elements(view)) {
      values[mapper.index(entity)] = f(entity);
    }

    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(view);
    vtkwriter.addCellData(values, "data");
    const std::string piecefilesFolderName = "piecefiles";
    const std::string piecefilesPath = f.dir() + "/" + piecefilesFolderName + "/";
    Common::test_create_directory(piecefilesPath);
    vtkwriter.pwrite(f.filename(), f.dir(), piecefilesFolderName, Dune::VTK::appendedraw);
  }

  template <class GridViewType>
  class FunctorBase
  {
  public:
    using Element = extract_entity_t<GridViewType>;
    FunctorBase(std::string filename = "Functor", const std::string dirname = ".")
      : filename_(filename)
      , dir_(dirname)
    {}

    virtual ~FunctorBase() {}

    const std::string filename() const
    {
      return filename_;
    }
    const std::string dir() const
    {
      return dir_;
    }

    virtual double operator()(const Element& /*ent*/) const = 0;

    std::vector<double> values(const GridViewType& view)
    {
      std::vector<double> ret(view.size(0));
      return ret;
    }

  protected:
    const std::string filename_;
    const std::string dir_;
  };

  template <class GridViewType>
  class VolumeFunctor : public FunctorBase<GridViewType>
  {
  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    VolumeFunctor(std::string filename = "VolumeFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
    {}

    double operator()(const Element& ent) const
    {
      return ent.geometry().volume();
    }
  };

  template <class GridViewType>
  class ProcessIdFunctor : public FunctorBase<GridViewType>
  {
  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    ProcessIdFunctor(std::string filename = "ProcessIDFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
    {}

    double operator()(const Element& /*ent*/) const
    {
      return Dune::MPIHelper::getCollectiveCommunication().rank();
    }
  };

  template <class GridViewType, bool enable = has_boundary_id<GridViewType>::value>
  class BoundaryIDFunctor : public FunctorBase<GridViewType>
  {
    const GridViewType& gridview_;

  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    BoundaryIDFunctor(const GridViewType& view,
                      std::string filename = "BoundaryIDFunctor",
                      const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
      , gridview_(view)
    {}

    double operator()(const Element& entity) const
    {
      double ret(0.0);
      int numberOfBoundarySegments(0);
      bool isOnBoundary = false;

      const auto intersection_it_end = gridview_.iend(entity);
      for (auto intersection_it = gridview_.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
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

  template <class GridViewType>
  class BoundaryIDFunctor<GridViewType, false> : public FunctorBase<GridViewType>
  {
    const GridViewType& gridview_;

  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    BoundaryIDFunctor(const GridViewType& view,
                      std::string filename = "BoundaryIDFunctor",
                      const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
      , gridview_(view)
    {
      DXTC_LOG_INFO_0 << "Boundary visualization for unsupported grid requested " << XT::Common::get_typename(gridview_)
                      << std::endl;
    }

    double operator()(const Element&) const
    {
      return -1;
    }
  };

  template <class GridViewType, class BoundaryInfoType>
  class BoundaryTypeFunctor : public FunctorBase<GridViewType>
  {
    const GridViewType& gridview_;
    const std::string type_;
    const BoundaryInfoType& boundaryInfo_;

  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    BoundaryTypeFunctor(const GridViewType& view,
                        const BoundaryInfoType& boundaryInfo,
                        std::string type,
                        std::string filename = "BoundaryTypeFunctor",
                        const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
      , gridview_(view)
      , type_(type)
      , boundaryInfo_(boundaryInfo)
    {}

    double operator()(const Element& entity) const
    {
      static const constexpr DirichletBoundary dirichlet_type{};
      static const constexpr NeumannBoundary neumann_type{};
      for (auto intersectionIt = gridview_.ibegin(entity); intersectionIt != gridview_.iend(entity); ++intersectionIt) {
        if (type_ == "dirichlet") {
          return (boundaryInfo_.type(*intersectionIt) == dirichlet_type);
        } else if (type_ == "neumann") {
          return (boundaryInfo_.type(*intersectionIt) == neumann_type);
        } else {
          DUNE_THROW(Common::Exceptions::internal_error, "Unknown type '" << type_ << "'!");
        }
      }
      return 0;
    }
  };

  template <class GridViewType>
  class AreaMarker : public FunctorBase<GridViewType>
  {

  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    AreaMarker(std::string filename = "AreaFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
    {}

    double operator()(const Element& entity) const
    {
      using EntityGeometryType = typename Element::Geometry;
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

  template <class GridViewType>
  class GeometryFunctor : public FunctorBase<GridViewType>
  {
  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    GeometryFunctor(std::string filename = "GeometryFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
    {}

    double operator()(const Element& ent) const
    {
      const auto& geo = ent.geometry();
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

  template <class GridViewType>
  class PartitionTypeFunctor : public FunctorBase<GridViewType>
  {
  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    PartitionTypeFunctor(std::string filename = "PartitionTypeFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
    {}

    double operator()(const Element& ent) const
    {
      const int type{static_cast<int>(ent.partitionType())};
      DXTC_LOG_ERROR << "TYPE " << type << std::endl;
      return static_cast<double>(type);
    }
  };

  template <class GridViewType, bool enable = has_boundary_id<GridViewType>::value>
  class IndexFunctor : public FunctorBase<GridViewType>
  {
    const GridViewType& gridview_;

  public:
    using Element = typename FunctorBase<GridViewType>::Element;
    IndexFunctor(const GridViewType& view, std::string filename = "IndexFunctor", const std::string dirname = ".")
      : FunctorBase<GridViewType>(filename, dirname)
      , gridview_(view)
    {}

    double operator()(const Element& entity) const
    {
      return gridview_.indexSet().index(entity);
    }
  };


  //! supply functor
  template <class Grid>
  static void all(const Grid& grid, const std::string outputDir = "visualisation")
  {
    // make function objects
    BoundaryIDFunctor<Grid> boundaryFunctor(grid, "boundaryFunctor", outputDir);
    AreaMarker<Grid> areaMarker("areaMarker", outputDir);
    GeometryFunctor<Grid> geometryFunctor("geometryFunctor", outputDir);
    ProcessIdFunctor<Grid> processIdFunctor("ProcessIdFunctor", outputDir);
    VolumeFunctor<Grid> volumeFunctor("volumeFunctor", outputDir);
    PartitionTypeFunctor<Grid> partitionTypeFunctor("partitionTypeFunctor", outputDir);

    // call the visualization functions
    const auto view = grid.leafGridView();
    elementdata(view, boundaryFunctor);
    elementdata(view, areaMarker);
    elementdata(view, geometryFunctor);
    elementdata(view, processIdFunctor);
    elementdata(view, volumeFunctor);
    elementdata(view, partitionTypeFunctor);
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

} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_OUTPUT_ENTITY_VISUALIZATION_HH
