// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr     (2014)
//   Barbara Verfürth (2015)
//   Felix Schindler  (2012 - 2016)
//   Kirsten Weber    (2012)
//   Rene Milk        (2012 - 2015)
//   Tobias Leibner   (2014)

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

#include <dune/xt/grid/structuredgridfactory.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/grid/grids.hh>

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

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID

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

#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID

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
  static GridProvider<GridType> create(const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
                                       const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
                                       const std::array<unsigned int, GridType::dimension>& num_elements,
                                       const unsigned int num_refinements,
                                       const std::array<unsigned int, GridType::dimension>& overlap_size)
  {
    static const int variant = ElementVariant<GridType>::id;
    static_assert(variant == 1 || variant == 2, "variant has to be 1 or 2!");
    for (unsigned int dd = 0; dd < GridType::dimension; ++dd) {
      if (!(lower_left[dd] < upper_right[dd]))
        DUNE_THROW(Common::Exceptions::wrong_input_given,
                   "lower_left has to be elementwise smaller than upper_right!\n\n"
                       << lower_left[dd]
                       << " vs. "
                       << upper_right[dd]);
    }
    std::shared_ptr<GridType> grd_ptr(nullptr);
    switch (variant) {
      case 1:
        grd_ptr = StructuredGridFactory<GridType>::createCubeGrid(lower_left, upper_right, num_elements, overlap_size);
        break;
      case 2:
      default:
        grd_ptr = StructuredGridFactory<GridType>::createSimplexGrid(lower_left, upper_right, num_elements);
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
    return GridProvider<GridType>(grd_ptr);
  } // ... create(...)

  static GridProvider<GridType>
  create(const typename GridType::ctype& lower_left,
         const typename GridType::ctype& upper_right,
         const unsigned int num_elements =
             cube_gridprovider_default_config().get<std::vector<unsigned int>>("num_elements")[0],
         const unsigned int num_refinements = cube_gridprovider_default_config().get<unsigned int>("num_refinements"),
         const unsigned int overlap_size =
             cube_gridprovider_default_config().get<std::vector<unsigned int>>("overlap_size")[0])
  {
    return create(FieldVector<typename GridType::ctype, GridType::dimension>(lower_left),
                  FieldVector<typename GridType::ctype, GridType::dimension>(upper_right),
                  Common::make_array<unsigned int, GridType::dimension>(num_elements),
                  num_refinements,
                  Common::make_array<unsigned int, GridType::dimension>(overlap_size));
  } // ... create(...)

  static GridProvider<GridType> create(const Common::Configuration& cfg = cube_gridprovider_default_config())
  {
    static const size_t d = GridType::dimension;
    auto overlap_size = cube_gridprovider_default_config().get<std::vector<unsigned int>>("overlap_size");
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
    auto lower_left =
        cfg.get("lower_left",
                cube_gridprovider_default_config().get<FieldVector<typename GridType::ctype, GridType::dimension>>(
                    "lower_left"));
    auto upper_right =
        cfg.get("upper_right",
                cube_gridprovider_default_config().get<FieldVector<typename GridType::ctype, GridType::dimension>>(
                    "upper_right"));
    auto num_elements =
        cfg.get("num_elements", cube_gridprovider_default_config().get<std::vector<unsigned int>>("num_elements"));
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
        cfg.get("num_refinements", cube_gridprovider_default_config().get<unsigned int>("num_refinements"));
    return create(lower_left, upper_right, num_elements_array, num_refinements, overlap_size_array);
  } // ... create(...)
}; // struct CubeGridProviderFactory


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType>>::type make_cube_grid(
    const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
    const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
    const std::array<unsigned int, GridType::dimension> num_elements =
        cube_gridprovider_default_config().get<std::vector<unsigned int>>("num_elements")[0],
    const unsigned int num_refinements = cube_gridprovider_default_config().get<unsigned int>("num_refinements"),
    const std::array<unsigned int, GridType::dimension> overlap_size =
        Common::make_array<unsigned int, GridType::dimension>(
            cube_gridprovider_default_config().get<std::vector<unsigned int>>("overlap_size")))
{
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size);
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType>>::type make_cube_grid(
    const typename GridType::ctype& lower_left,
    const typename GridType::ctype& upper_right,
    const unsigned int num_elements =
        cube_gridprovider_default_config().get<std::vector<unsigned int>>("num_elements")[0],
    const unsigned int num_refinements = cube_gridprovider_default_config().get<unsigned int>("num_refinements"),
    const unsigned int overlap_size =
        cube_gridprovider_default_config().get<std::vector<unsigned int>>("overlap_size")[0])
{
  return CubeGridProviderFactory<GridType>::create(
      lower_left, upper_right, num_elements, num_refinements, overlap_size);
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType>>::type
make_cube_grid(const Common::Configuration& cfg = cube_gridprovider_default_config())
{
  return CubeGridProviderFactory<GridType>::create(cfg);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_CUBE_HH
