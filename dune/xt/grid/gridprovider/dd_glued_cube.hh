// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr     (2014)
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2017, 2019)
//   Kirsten Weber    (2012)
//   René Fritze      (2012 - 2019)
//   Tobias Leibner   (2014, 2016, 2018, 2020)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_DD_GLUED_CUBE_HH
#define DUNE_XT_GRID_GRIDPROVIDER_DD_GLUED_CUBE_HH

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
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/structuredgridfactory.hh>
#include <dune/xt/grid/type_traits.hh>

#include "provider.hh"
#include "dune/xt/grid/dd/glued.hh"

namespace Dune::XT::Grid {


static inline std::string cube_glued_id()
{
  return "xt.grid.dd.glued";
}


template <class MacroGridType, class LocalGridType>
class CubeGluedGridProviderFactory
{
//  static_assert(is_grid<GluedGridType>::value);

  template <typename G>
  struct ElementVariant;

  template <typename G>
  struct ElementVariant
  {
    static constexpr int id = 2;
  };

  template <int dim, class Coords>
  struct ElementVariant<Dune::YaspGrid<dim, Coords>>
  {
    static constexpr int id = 1;
  };

#if HAVE_DUNE_SPGRID
  template <class ct, int dim, template <int> class Refinement, class Comm>
  struct ElementVariant<Dune::SPGrid<ct, dim, Refinement, Comm>>
  {
    static constexpr int id = 1;
  };
#endif

#if HAVE_DUNE_ALUGRID

  template <int dimGrid, int dimWorld, class MpiCommImp>
  struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::conforming, MpiCommImp>>
  {
    static constexpr int id = 1;
  };

  template <int dimGrid, int dimWorld, class MpiCommImp>
  struct ElementVariant<Dune::ALUGrid<dimGrid, dimWorld, Dune::cube, Dune::nonconforming, MpiCommImp>>
  {
    static constexpr int id = 1;
  };

#endif // HAVE_DUNE_ALUGRID

public:
  static constexpr bool available = true;

  static std::string static_id()
  {
    return cube_glued_id();
  }

  /// TODO simplex grid overlap_size
  static XT::Grid::DD::Glued<MacroGridType, LocalGridType, XT::Grid::Layers::leaf> create(XT::Grid::GridProvider<MacroGridType>& macro_grid,
                                                                                          const unsigned int num_refinements)
  {
//    static constexpr int variant = ElementVariant<MacroGridType>::id;
//    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
//    for (unsigned int dd = 0; dd < MacroGridType::dimension; ++dd) {
//      if (!(lower_left[dd] < upper_right[dd]))
//        DUNE_THROW(Common::Exceptions::wrong_input_given,
//                   "lower_left has to be elementwise smaller than upper_right!\n\nlower_left = "
//                       << lower_left << "\n\nupper_right = " << upper_right);
//    }
    std::shared_ptr<MacroGridType> grd_ptr(nullptr);
    switch (variant) {
      case 1:
        grd_ptr = XT::Grid::StructuredGridFactory<MacroGridType>::createCubeGrid(
            lower_left, upper_right, num_elements, overlap_size, mpi_comm);
        break;
      case 2:
      default:
        grd_ptr = XT::Grid::StructuredGridFactory<MacroGridType>::createSimplexGrid(
            lower_left, upper_right, num_elements, mpi_comm);
        break;
    }
    grd_ptr->loadBalance();
#if HAVE_ALBERTA
    if (!std::is_same<MacroGridType, AlbertaGrid<MacroGridType::dimension, MacroGridType::dimension>>::value)
#endif
      grd_ptr->preAdapt();
    grd_ptr->globalRefine(boost::numeric_cast<int>(num_refinements));
#if HAVE_ALBERTA
    if (!std::is_same<MacroGridType, AlbertaGrid<MacroGridType::dimension, MacroGridType::dimension>>::value)
#endif
      grd_ptr->postAdapt();
    grd_ptr->loadBalance();
    return GridProvider<MacroGridType>(grd_ptr);
  } // ... create(...)

  static GridProvider<MacroGridType> create(const typename MacroGridType::ctype& lower_left,
                                       const typename MacroGridType::ctype& upper_right,
                                       const unsigned int num_elements,
                                       const unsigned int num_refinements,
                                       const unsigned int overlap_size,
                                       MPIHelper::MPICommunicator mpi_comm)
  {
    return create(FieldVector<typename MacroGridType::ctype, MacroGridType::dimension>(lower_left),
                  FieldVector<typename MacroGridType::ctype, MacroGridType::dimension>(upper_right),
                  Common::make_array<unsigned int, MacroGridType::dimension>(num_elements),
                  num_refinements,
                  Common::make_array<unsigned int, MacroGridType::dimension>(overlap_size),
                  mpi_comm);
  } // ... create(...)

  static GridProvider<MacroGridType> create(const Common::Configuration& cfg, MPIHelper::MPICommunicator mpi_comm)
  {
    static constexpr size_t d = MacroGridType::dimension;
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
                 "overlap_size has to be a single number or a vector with at least "
                     << d << "elements, has only " << overlap_size.size() << " elements!");
    }
    auto lower_left = cfg.get(
        "lower_left",
        cube_gridprovider_default_config().template get<FieldVector<typename MacroGridType::ctype, MacroGridType::dimension>>(
            "lower_left"));
    auto upper_right = cfg.get(
        "upper_right",
        cube_gridprovider_default_config().template get<FieldVector<typename MacroGridType::ctype, MacroGridType::dimension>>(
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
                 "num_elements has to be a single number or a vector with at least "
                     << d << "elements, has only " << num_elements.size() << " elements!");
    }
    auto num_refinements =
        cfg.get("num_refinements", cube_gridprovider_default_config().template get<unsigned int>("num_refinements"));
    return create(lower_left, upper_right, num_elements_array, num_refinements, overlap_size_array, mpi_comm);
  } // ... create(...)
}; // struct CubeGluedGridProviderFactory


template <class GridType>
auto make_cube_grid(
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
  static_assert(is_grid<GridType>::value);
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm);
}


template <class GridType>
auto make_cube_grid(
    const typename GridType::ctype& lower_left,
    const typename GridType::ctype& upper_right,
    const unsigned int num_elements =
        cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements").at(0),
    const unsigned int num_refinements =
        cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
    const unsigned int overlap_size =
        cube_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size").at(0),
    MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm);
}


template <class GridType>
auto make_cube_grid(const Common::Configuration& cfg = cube_gridprovider_default_config(),
                    MPIHelper::MPICommunicator mpi_comm = MPIHelper::getCommunicator())
{
  static_assert(is_grid<GridType>::value);
  return CubeGridProviderFactory<GridType>::create(cfg, mpi_comm);
}


} // namespace Dune::XT::Grid


#endif // DUNE_XT_GRID_GRIDPROVIDER_CUBE_HH
