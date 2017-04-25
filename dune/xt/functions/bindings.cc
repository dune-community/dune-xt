// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/xt/common/exceptions.hh>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include <dune/xt/grid/grids.hh>

#include "interfaces.pbh"
#include "constant.pbh"
#include "checkerboard.pbh"
#include "expression.pbh"
#include "spe10.pbh"


template <class G>
void addbind_for_Grid(pybind11::module& m, const std::string& grid_id)
{
  auto i_1_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 1, 1>(m, grid_id);
  auto i_2_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 2, 1>(m, grid_id);
  auto i_3_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 3, 1>(m, grid_id);
  auto i_4_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 4, 1>(m, grid_id);

  auto i_2_2 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 2, 2>(m, grid_id);

  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         1,
                                                         1,
                                                         1,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         2,
                                                         1,
                                                         2,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       2,
                                                       1,
                                                       2,
                                                       1>(i_2_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         3,
                                                         1,
                                                         3,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       3,
                                                       1,
                                                       3,
                                                       1>(i_3_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         4,
                                                         1,
                                                         4,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       4,
                                                       1,
                                                       4,
                                                       1>(i_4_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         2,
                                                         2,
                                                         2,
                                                         2>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       2,
                                                       2,
                                                       2,
                                                       2>(i_2_2);


  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, G::dimension, Dune::XT::Functions::internal::Combination::sum, 1, 1, 1, 1>(
          m, grid_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<G,
                                                                        G::dimension,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        1,
                                                                        1,
                                                                        1,
                                                                        1>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, G::dimension, Dune::XT::Functions::internal::Combination::sum, 2, 1, 2, 1>(
          m, grid_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<G,
                                                                        G::dimension,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        2,
                                                                        1,
                                                                        2,
                                                                        1>(i_2_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, G::dimension, Dune::XT::Functions::internal::Combination::sum, 3, 1, 3, 1>(
          m, grid_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<G,
                                                                        G::dimension,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        3,
                                                                        1,
                                                                        3,
                                                                        1>(i_3_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, G::dimension, Dune::XT::Functions::internal::Combination::sum, 4, 1, 4, 1>(
          m, grid_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<G,
                                                                        G::dimension,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        4,
                                                                        1,
                                                                        4,
                                                                        1>(i_4_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, G::dimension, Dune::XT::Functions::internal::Combination::sum, 2, 2, 2, 2>(
          m, grid_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<G,
                                                                        G::dimension,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        2,
                                                                        2,
                                                                        2,
                                                                        2>(i_2_2);

  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         1,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         2,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       2,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         3,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       3,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         4,
                                                         1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       4,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<G,
                                                         G::dimension,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         2,
                                                         2>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       G::dimension,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       2,
                                                       2>(i_1_1);

  Dune::XT::Functions::bind_ConstantFunction<G, G::dimension, 1, 1>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, G::dimension, 2, 1>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, G::dimension, 3, 1>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, G::dimension, 4, 1>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, G::dimension, 2, 2>(m, grid_id);

  Dune::XT::Functions::bind_CheckerboardFunction<G, G::dimension, 1, 1>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, G::dimension, 2, 1>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, G::dimension, 3, 1>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, G::dimension, 4, 1>(m, grid_id);

  Dune::XT::Functions::bind_ExpressionFunction<G, G::dimension, 1, 1>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, G::dimension, 2, 1>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, G::dimension, 3, 1>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, G::dimension, 4, 1>(m, grid_id);

  Dune::XT::Functions::bind_Spe10Model1Function<G, G::dimension, 1, 1>(m, grid_id);
} // ... addbind_for_Grid(...)


PYBIND11_PLUGIN(_functions)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("_functions", "dune-xt-functions");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");

  addbind_for_Grid<Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>>(m, "1d_yaspgrid");
  addbind_for_Grid<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>(m, "2d_cube_yaspgrid");
#if HAVE_DUNE_ALUGRID
  addbind_for_Grid<Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>>(m, "2d_simplex_aluconform");
#endif
  //#if HAVE_UG
  //  addbind_for_Grid<Dune::UGGrid<2>>(m, "2d_simplex_uggrid");
  //#endif
  //#if HAVE_ALBERTA
  //  addbind_for_Grid<Dune::AlbertaGrid<2, 2>>(m, "2d_simplex_albertagrid");
  //#endif

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = boost::numeric_cast<int>(args.size());
          char** argv = Dune::XT::Common::vector_to_main_args(args);
          Dune::MPIHelper::instance(argc, argv);
#if HAVE_DUNE_FEM
          Dune::Fem::MPIManager::initialize(argc, argv);
#endif
        },
        "args"_a = std::vector<std::string>());

  m.def("_init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = -1,
        "max_debug_level"_a = -1,
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  return m.ptr();
}


#endif // HAVE_DUNE_PYBINDXI
