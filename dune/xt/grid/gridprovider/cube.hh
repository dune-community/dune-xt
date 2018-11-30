// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr     (2014)
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2017)
//   Kirsten Weber    (2012)
//   René Fritze      (2012 - 2018)
//   Tobias Leibner   (2014, 2016, 2018)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_CUBE_HH
#define DUNE_XT_GRID_GRIDPROVIDER_CUBE_HH

#include <array>
#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>
#include <limits>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/misc.hh>
#include <dune/xt/grid/dd/subdomains/factory.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/structuredgridfactory.hh>
#include <dune/xt/grid/type_traits.hh>

#include "provider.hh"

namespace Dune {
namespace XT {
namespace Grid {


static inline std::string cube_gridprovider_id()
{
  return "xt.grid.gridprovider.cube";
}


static inline Common::Configuration cube_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"] = cube_gridprovider_id();
  config["lower_left"] = "[0 0 0 0]";
  config["upper_right"] = "[1 1 1 1]";
  config["num_elements"] = "[8 8 8 8]";
  config["num_refinements"] = "0";
  config["overlap_size"] = "[1 1 1 1]";
  return config;
}


static inline std::string cube_dd_subdomains_gridprovider_id()
{
  return "xt.grid.gridprovider.dd.subdomains.cube";
}


static inline Common::Configuration cube_dd_subdomains_gridprovider_default_config()
{
  Common::Configuration config = cube_gridprovider_default_config();
  config["num_partitions"] = "[2 2 2 2]";
  config["oversampling_layers"] = "0";
  config["inner_boundary_segment_index"] = Common::to_string(std::numeric_limits<size_t>::max() - 42);
  return config;
}


template <class GridType>
class CubeGridProviderFactory
{
  static_assert(is_grid<GridType>::value, "");

  template <typename G>
  struct ElementVariant;

  template <typename G>
  struct ElementVariant
  {
    static const int id = 2;
  };

  template <int dim, class Coords>
  struct ElementVariant<Dune::YaspGrid<dim, Coords>>
  {
    static const int id = 1;
  };

#if HAVE_DUNE_SPGRID
  template <class ct, int dim, template <int> class Refinement, class Comm>
  struct ElementVariant<Dune::SPGrid<ct, dim, Refinement, Comm>>
  {
    static const int id = 1;
  };
#endif

#if HAVE_DUNE_ALUGRID

  template <int dimGrid, int dimWorld, class MpiCommImp>
  struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::conforming, MpiCommImp>>
  {
    static const int id = 1;
  };

  template <int dimGrid, int dimWorld, class MpiCommImp>
  struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::nonconforming, MpiCommImp>>
  {
    static const int id = 1;
  };

#endif // HAVE_DUNE_ALUGRID

public:
  static const bool available = true;

  static std::string static_id()
  {
    return cube_gridprovider_id();
  }

  static Common::Configuration default_config()
  {
    return cube_gridprovider_default_config();
  }

  /// TODO simplex grid overlap_size
  static GridProvider<GridType, none_t>
  create(const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
         const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
         const std::array<unsigned int, GridType::dimension>& num_elements,
         const unsigned int num_refinements,
         const std::array<unsigned int, GridType::dimension>& overlap_size,
         MPIHelper::MPICommunicator mpi_comm)
  {
    static const int variant = ElementVariant<GridType>::id;
    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
    for (unsigned int dd = 0; dd < GridType::dimension; ++dd) {
      if (!(lower_left[dd] < upper_right[dd]))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "lower_left has to be elementwise smaller than upper_right!\n\nlower_left = " << lower_left
                                                                                                 << "\n\nupper_right = "
                                                                                                 << upper_right);
    }
    std::shared_ptr<GridType> grd_ptr(nullptr);
    switch (variant) {
      case 1:
        grd_ptr = XT::Grid::StructuredGridFactory<GridType>::createCubeGrid(
            lower_left, upper_right, num_elements, overlap_size, mpi_comm);
        break;
      case 2:
      default:
        grd_ptr = XT::Grid::StructuredGridFactory<GridType>::createSimplexGrid(
            lower_left, upper_right, num_elements, mpi_comm);
        break;
    }
    grd_ptr->loadBalance();
#if HAVE_ALBERTA
    if (!std::is_same<GridType, AlbertaGrid<GridType::dimension, GridType::dimension>>::value)
#endif
      grd_ptr->preAdapt();
    grd_ptr->globalRefine(boost::numeric_cast<int>(num_refinements));
#if HAVE_ALBERTA
    if (!std::is_same<GridType, AlbertaGrid<GridType::dimension, GridType::dimension>>::value)
#endif
      grd_ptr->postAdapt();
    grd_ptr->loadBalance();
    return GridProvider<GridType, none_t>(grd_ptr);
  } // ... create(...)

  static GridProvider<GridType, none_t> create(const typename GridType::ctype& lower_left,
                                               const typename GridType::ctype& upper_right,
                                               const unsigned int num_elements,
                                               const unsigned int num_refinements,
                                               const unsigned int overlap_size,
                                               MPIHelper::MPICommunicator mpi_comm)
  {
    return create(FieldVector<typename GridType::ctype, GridType::dimension>(lower_left),
                  FieldVector<typename GridType::ctype, GridType::dimension>(upper_right),
                  Common::make_array<unsigned int, GridType::dimension>(num_elements),
                  num_refinements,
                  Common::make_array<unsigned int, GridType::dimension>(overlap_size),
                  mpi_comm);
  } // ... create(...)

  static GridProvider<GridType, none_t> create(const Common::Configuration& cfg, MPIHelper::MPICommunicator mpi_comm)
  {
    static const size_t d = GridType::dimension;
    auto overlap_size =
        cfg.has_key("overlap_size")
            ? cfg.template get<std::vector<unsigned int>>("overlap_size")
            : cube_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size");
    std::array<unsigned int, d> overlap_size_array;
    if (overlap_size.size() >= d) {
      overlap_size_array = Common::make_array<unsigned int, d>(overlap_size);
    } else if (overlap_size.size() >= 1) {
      overlap_size_array = Common::make_array<unsigned int, d>(overlap_size);
    } else {
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "overlap_size has to be a single number or a vector with at least " << d << "elements, has only "
                                                                                     << overlap_size.size()
                                                                                     << " elements!");
    }
    auto lower_left = cfg.get(
        "lower_left",
        cube_gridprovider_default_config().template get<FieldVector<typename GridType::ctype, GridType::dimension>>(
            "lower_left"));
    auto upper_right = cfg.get(
        "upper_right",
        cube_gridprovider_default_config().template get<FieldVector<typename GridType::ctype, GridType::dimension>>(
            "upper_right"));
    auto num_elements = cfg.get(
        "num_elements", cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements"));
    std::array<unsigned int, d> num_elements_array;
    if (num_elements.size() >= d) {
      num_elements_array = Common::make_array<unsigned int, d>(num_elements);
    } else if (num_elements.size() >= 1) {
      num_elements_array = Common::make_array<unsigned int, d>(num_elements);
    } else {
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "num_elements has to be a single number or a vector with at least " << d << "elements, has only "
                                                                                     << num_elements.size()
                                                                                     << " elements!");
    }
    auto num_refinements =
        cfg.get("num_refinements", cube_gridprovider_default_config().template get<unsigned int>("num_refinements"));
    return create(lower_left, upper_right, num_elements_array, num_refinements, overlap_size_array, mpi_comm);
  } // ... create(...)
}; // struct CubeGridProviderFactory


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType, none_t>>::type make_cube_grid(
    const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
    const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
    const std::array<unsigned int, GridType::dimension> num_elements =
        cube_gridprovider_default_config().template get<std::array<unsigned int, GridType::dimension>>("num_elements"),
    const unsigned int num_refinements =
        cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
    const std::array<unsigned int, GridType::dimension> overlap_size =
        cube_gridprovider_default_config().template get<std::array<unsigned int, GridType::dimension>>("overlap_size"),
    MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm);
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType, none_t>>::type
make_cube_grid(const typename GridType::ctype& lower_left,
               const typename GridType::ctype& upper_right,
               const unsigned int num_elements =
                   cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements").at(0),
               const unsigned int num_refinements =
                   cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
               const unsigned int overlap_size =
                   cube_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size").at(0),
               MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm);
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType, none_t>>::type
make_cube_grid(const Common::Configuration& cfg = cube_gridprovider_default_config(),
               MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  return CubeGridProviderFactory<GridType>::create(cfg, mpi_comm);
}


template <class GridType>
class CubeDdSubdomainsGridProviderFactory
{
  static_assert(XT::Grid::is_grid<GridType>::value, "");

public:
  typedef DD::SubdomainGrid<GridType> DdGridType;

  static const bool available = true;

  static std::string static_id()
  {
    return cube_dd_subdomains_gridprovider_id();
  }

  static XT::Common::Configuration default_config()
  {
    return cube_dd_subdomains_gridprovider_default_config();
  }

  static GridProvider<GridType, DdGridType> create(const Common::Configuration& cfg,
                                                   MPIHelper::MPICommunicator mpi_comm)
  {
    const auto default_cfg = cube_dd_subdomains_gridprovider_default_config();
    return create(
        cfg.get("lower_left",
                default_cfg.get<FieldVector<typename GridType::ctype, GridType::dimension>>("lower_left")),
        cfg.get("upper_right",
                default_cfg.get<FieldVector<typename GridType::ctype, GridType::dimension>>("upper_right")),
        cfg.get("num_elements", default_cfg.get<std::array<unsigned int, GridType::dimension>>("num_elements")),
        cfg.get("num_refinements", default_cfg.get<unsigned int>("num_refinements")),
        cfg.get("overlap_size", default_cfg.get<std::array<unsigned int, GridType::dimension>>("overlap_size")),
        cfg.get("num_partitions", default_cfg.get<std::array<unsigned int, GridType::dimension>>("num_partitions")),
        cfg.get("oversampling_layers", default_cfg.get<size_t>("oversampling_layers")),
        cfg.get("inner_boundary_segment_index", default_cfg.get<size_t>("inner_boundary_segment_index")),
        mpi_comm);
  }

  static GridProvider<GridType, DdGridType>
  create(const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
         const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
         const std::array<unsigned int, GridType::dimension>& num_elements,
         const unsigned int num_refinements,
         const std::array<unsigned int, GridType::dimension>& overlap_size,
         const std::array<unsigned int, GridType::dimension>& num_partitions,
         const size_t num_oversampling_layers,
         const size_t inner_boundary_segment_index,
         MPIHelper::MPICommunicator mpi_comm)
  {
    auto grid = make_cube_grid<GridType>(lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm)
                    .grid_ptr();
    grid->loadBalance();
    const auto dims = XT::Grid::dimensions(grid->leafGridView());
    const auto bounding_box = dims.bounding_box();

    typedef DD::SubdomainGridFactory<GridType> DdGridFactoryType;
    const size_t neighbor_recursion_level = DD::internal::NeighborRecursionLevel<GridType>::compute();
    // prepare
    DdGridFactoryType factory(*grid, inner_boundary_segment_index);
    factory.prepare();
    // global grid part
    const auto global_grid_part = factory.globalGridView();
    // walk the grid
    const auto entity_it_end = global_grid_part->template end<0, All_Partition>();
    for (auto entity_it = global_grid_part->template begin<0, All_Partition>(); entity_it != entity_it_end;
         ++entity_it) {
      // get center of entity
      const auto& entity = *entity_it;
      const auto center = entity.geometry().center();
      // decide on the subdomain this entity shall belong to
      std::vector<size_t> whichPartition(GridType::dimension, 0);
      for (size_t dd = 0; dd < GridType::dimension; ++dd)
        whichPartition[dd] =
            (std::min((unsigned int)(std::floor(num_partitions[dd] * ((center[dd] - bounding_box[0][dd])
                                                                      / (bounding_box[1][dd] - bounding_box[0][dd])))),
                      num_partitions[dd] - 1));
      size_t subdomain = 0;
      if (GridType::dimension == 1)
        subdomain = whichPartition[0];
      else if (GridType::dimension == 2)
        subdomain = whichPartition[0] + whichPartition[1] * num_partitions[0];
      else if (GridType::dimension == 3)
        subdomain = whichPartition[0] + whichPartition[1] * num_partitions[0]
                    + whichPartition[2] * num_partitions[1] * num_partitions[0];
      else
        DUNE_THROW(Dune::NotImplemented,
                   "ERROR in " << static_id() << ": not implemented for grid dimensions other than 1, 2 or 3!");
      // add entity to subdomain
      factory.add(entity, subdomain /*, prefix + "  ", out*/);
    } // walk the grid
    // finalize
    factory.finalize(num_oversampling_layers, neighbor_recursion_level /*, prefix + "  ", out*/);
    // be done with it
    return GridProvider<GridType, DdGridType>(grid, factory.createMsGrid());
  } // ... create(...)
}; // class CubeDdSubdomainsGridProviderFactory


template <class GridType>
typename std::enable_if<XT::Grid::is_grid<GridType>::value, GridProvider<GridType, DD::SubdomainGrid<GridType>>>::type
make_cube_dd_subdomains_grid(
    const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
    const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
    const std::array<unsigned int, GridType::dimension> num_elements =
        cube_dd_subdomains_gridprovider_default_config().template get<std::array<unsigned int, GridType::dimension>>(
            "num_elements"),
    const unsigned int num_refinements =
        cube_dd_subdomains_gridprovider_default_config().template get<unsigned int>("num_refinements"),
    const std::array<unsigned int, GridType::dimension> overlap_size =
        cube_dd_subdomains_gridprovider_default_config().template get<std::array<unsigned int, GridType::dimension>>(
            "overlap_size"),
    const std::array<unsigned int, GridType::dimension> num_partitions =
        cube_dd_subdomains_gridprovider_default_config().template get<std::array<unsigned int, GridType::dimension>>(
            "num_partitions"),
    const size_t num_oversampling_layers =
        cube_dd_subdomains_gridprovider_default_config().template get<size_t>("num_refinements"),
    const size_t inner_boundary_segment_index =
        cube_dd_subdomains_gridprovider_default_config().template get<size_t>("inner_boundary_segment_index"),
    MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  return CubeDdSubdomainsGridProviderFactory<GridType>::create(lower_left,
                                                               upper_right,
                                                               num_elements,
                                                               num_refinements,
                                                               overlap_size,
                                                               num_partitions,
                                                               num_oversampling_layers,
                                                               inner_boundary_segment_index,
                                                               mpi_comm);
}

template <class GridType>
typename std::enable_if<XT::Grid::is_grid<GridType>::value, GridProvider<GridType, DD::SubdomainGrid<GridType>>>::type
make_cube_dd_subdomains_grid(const Common::Configuration& cfg = cube_dd_subdomains_gridprovider_default_config(),
                             MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  return CubeDdSubdomainsGridProviderFactory<GridType>::create(cfg, mpi_comm);
}


} // namespace Grid
} // namespace XT
} // namespace Dune


#include "cube.lib.hh"


#endif // DUNE_XT_GRID_GRIDPROVIDER_CUBE_HH
