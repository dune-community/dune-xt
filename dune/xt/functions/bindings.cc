#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>


#include <dune/xt/grid/grids.hh>

#include "interfaces.pbh"
#include "constant.pbh"
#include "checkerboard.pbh"
#include "expression.pbh"
#include "spe10.pbh"

namespace py = pybind11;
using namespace pybind11::literals;


PYBIND11_PLUGIN(functions)
{
  py::module m("functions", "dune-xt-functions");

  py::module::import("common");
  py::module::import("grid");

  typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> G;
  const std::string grid_id = "2d_simplex_aluconform";

  auto i_1_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 1, 1>(m, grid_id);
  auto i_1_2 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 1, 2>(m, grid_id);
  auto i_1_3 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 1, 3>(m, grid_id);
  auto i_1_4 = Dune::XT::Functions::bind_LocalizableFunctionInterface<G, 1, 4>(m, grid_id);

  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::difference, 1, 1, 1, 1>(m,
                                                                                                               grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::difference, 1, 2, 1, 2>(m,
                                                                                                               grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       2,
                                                       1,
                                                       2>(i_1_2);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::difference, 1, 3, 1, 3>(m,
                                                                                                               grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       3,
                                                       1,
                                                       3>(i_1_3);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::difference, 1, 4, 1, 4>(m,
                                                                                                               grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       4,
                                                       1,
                                                       4>(i_1_4);

  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::sum, 1, 1, 1, 1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G, Dune::XT::Functions::internal::Combination::sum, 1, 1, 1, 1>(
          i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::sum, 1, 2, 1, 2>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G, Dune::XT::Functions::internal::Combination::sum, 1, 2, 1, 2>(
          i_1_2);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::sum, 1, 3, 1, 3>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G, Dune::XT::Functions::internal::Combination::sum, 1, 3, 1, 3>(
          i_1_3);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::sum, 1, 4, 1, 4>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G, Dune::XT::Functions::internal::Combination::sum, 1, 4, 1, 4>(
          i_1_4);

  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::product, 1, 1, 1, 1>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::product, 1, 1, 1, 2>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       2>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::product, 1, 1, 1, 3>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       3>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<G, Dune::XT::Functions::internal::Combination::product, 1, 1, 1, 4>(m, grid_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<G,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       4>(i_1_1);

  Dune::XT::Functions::bind_ConstantFunction<G, 1, 1>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, 1, 2>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, 1, 3>(m, grid_id);
  Dune::XT::Functions::bind_ConstantFunction<G, 1, 4>(m, grid_id);

  Dune::XT::Functions::bind_CheckerboardFunction<G, 1, 1>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, 1, 2>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, 1, 3>(m, grid_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<G, 1, 4>(m, grid_id);

  Dune::XT::Functions::bind_ExpressionFunction<G, 1, 1>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, 1, 2>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, 1, 3>(m, grid_id);
  Dune::XT::Functions::bind_ExpressionFunction<G, 1, 4>(m, grid_id);

  Dune::XT::Functions::bind_Spe10Model1Function<G, 1, 1>(m, grid_id);

  return m.ptr();
}


#endif // HAVE_DUNE_PYBINDXI
