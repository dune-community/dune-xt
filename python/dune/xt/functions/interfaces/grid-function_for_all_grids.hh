// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_FOR_ALL_GRIDS_HH
#define PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_FOR_ALL_GRIDS_HH

#include <dune/xt/grid/type_traits.hh>

#include <python/dune/xt/functions/base/combined-grid-function.hh>

#include "grid-function.hh"


template <size_t ii>
struct Int
{
  static const constexpr size_t value = ii;
};


template <class GridTypes>
struct GridFunctionInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;

  template <size_t r, class Dims = boost::tuple<Int<1>, Int<2>, Int<3>>>
  struct for_all_rC
  {
    static const constexpr size_t rC = Dims::head_type::value;

    template <bool vector = (r != 1 && rC == 1), bool matrix = (rC != 1), bool anything = false>
    struct product_helper // <true, false, ...>
    {
      static void addbind(pybind11::module& m)
      {
        using Dune::XT::Functions::bindings::ProductGridFunction;
        using Dune::XT::Grid::bindings::grid_name;

        // special case: vector * vector
        ProductGridFunction<G, E, r, 1, r, 1>::bind(m, grid_name<G>::value());
      }
    };

    template <bool anything>
    struct product_helper<false, true, anything>
    {
      static void addbind(pybind11::module& m)
      {
        using Dune::XT::Functions::bindings::ProductGridFunction;
        using Dune::XT::Grid::bindings::grid_name;

        // general case: matrix * matrix or vector
        ProductGridFunction<G, E, r, rC, rC, 1>::bind(m, grid_name<G>::value());
        ProductGridFunction<G, E, r, rC, rC, 2>::bind(m, grid_name<G>::value());
        ProductGridFunction<G, E, r, rC, rC, 3>::bind(m, grid_name<G>::value());
      }
    };

    template <bool scalar = (r == 1 && rC == 1), bool anything = true>
    struct fraction_helper
    {
      static void addbind(pybind11::module& m)
      {
        using Dune::XT::Functions::bindings::FractionGridFunction;
        using Dune::XT::Grid::bindings::grid_name;

        FractionGridFunction<G, E>::bind(m, grid_name<G>::value());
      }
    };

    template <bool a>
    struct fraction_helper<false, a>
    {
      static void addbind(pybind11::module& /*m*/) {}
    };

    static void bind(pybind11::module& m)
    {
      using Dune::XT::Functions::bindings::DifferenceGridFunction;
      using Dune::XT::Functions::bindings::GridFunctionInterface;
      using Dune::XT::Functions::bindings::ProductGridFunction;
      using Dune::XT::Functions::bindings::SumGridFunction;
      using Dune::XT::Grid::bindings::grid_name;

      GridFunctionInterface<G, E, r, rC>::bind(m, grid_name<G>::value());
      DifferenceGridFunction<G, E, r, rC>::bind(m, grid_name<G>::value());
      SumGridFunction<G, E, r, rC>::bind(m, grid_name<G>::value());
      // we can always multiply with a scalar from the right ...
      ProductGridFunction<G, E, r, rC, 1, 1>::bind(m, grid_name<G>::value());
      // .. and with lots of other dims
      product_helper<>::addbind(m);
      fraction_helper<>::addbind(m);

      for_all_rC<r, typename Dims::tail_type>::bind(m);
    }
  };

  template <size_t r>
  struct for_all_rC<r, boost::tuples::null_type>
  {
    static void bind(pybind11::module& /*m*/) {}
  };


  template <class Dims = boost::tuple<Int<1>, Int<2>, Int<3>>, bool anything = false>
  struct for_all_r_and_rC
  {
    static const constexpr size_t r = Dims::head_type::value;

    static void bind(pybind11::module& m)
    {
      for_all_rC<r>::bind(m);

      for_all_r_and_rC<typename Dims::tail_type>::bind(m);
    }
  };

  template <bool a>
  struct for_all_r_and_rC<boost::tuples::null_type, a>
  {
    static void bind(pybind11::module& /*m*/) {}
  };


  static void bind(pybind11::module& m)
  {
    using Dune::XT::Functions::bindings::FractionGridFunction;
    using Dune::XT::Grid::bindings::grid_name;

    for_all_r_and_rC<>::bind(m);

    GridFunctionInterface_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct GridFunctionInterface_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


#endif // PYTHON_DUNE_XT_FUNCTIONS_INTERFACES_GRID_FUNCTION_FOR_ALL_GRIDS_HH
