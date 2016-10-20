#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

// PYBIND11_MAKE_OPAQUE(std::vector<ssize_t>);
// PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
// PYBIND11_MAKE_OPAQUE(std::vector<double>);

// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<ssize_t>>);
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::string>>);
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);

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

  typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> YaspGrid2dType;
  const std::string yasp_grid_2d_id = "2d_cube_yaspgrid";

  auto i_1_1 = Dune::XT::Functions::bind_LocalizableFunctionInterface<YaspGrid2dType, 1, 1>(m, yasp_grid_2d_id);
  auto i_1_2 = Dune::XT::Functions::bind_LocalizableFunctionInterface<YaspGrid2dType, 1, 2>(m, yasp_grid_2d_id);
  auto i_1_3 = Dune::XT::Functions::bind_LocalizableFunctionInterface<YaspGrid2dType, 1, 3>(m, yasp_grid_2d_id);
  auto i_1_4 = Dune::XT::Functions::bind_LocalizableFunctionInterface<YaspGrid2dType, 1, 4>(m, yasp_grid_2d_id);

  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         1,
                                                         1,
                                                         1,
                                                         1>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         1,
                                                         2,
                                                         1,
                                                         2>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       2,
                                                       1,
                                                       2>(i_1_2);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         1,
                                                         3,
                                                         1,
                                                         3>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       3,
                                                       1,
                                                       3>(i_1_3);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::difference,
                                                         1,
                                                         4,
                                                         1,
                                                         4>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::difference,
                                                       1,
                                                       4,
                                                       1,
                                                       4>(i_1_4);

  Dune::XT::Functions::
      bind_combined_LocalizableFunction<YaspGrid2dType, Dune::XT::Functions::internal::Combination::sum, 1, 1, 1, 1>(
          m, yasp_grid_2d_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        1,
                                                                        1,
                                                                        1,
                                                                        1>(i_1_1);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<YaspGrid2dType, Dune::XT::Functions::internal::Combination::sum, 1, 2, 1, 2>(
          m, yasp_grid_2d_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        1,
                                                                        2,
                                                                        1,
                                                                        2>(i_1_2);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<YaspGrid2dType, Dune::XT::Functions::internal::Combination::sum, 1, 3, 1, 3>(
          m, yasp_grid_2d_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        1,
                                                                        3,
                                                                        1,
                                                                        3>(i_1_3);
  Dune::XT::Functions::
      bind_combined_LocalizableFunction<YaspGrid2dType, Dune::XT::Functions::internal::Combination::sum, 1, 4, 1, 4>(
          m, yasp_grid_2d_id);
  Dune::XT::Functions::addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                                        Dune::XT::Functions::internal::Combination::sum,
                                                                        1,
                                                                        4,
                                                                        1,
                                                                        4>(i_1_4);

  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         1,
                                                         1>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       1>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         1,
                                                         2>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       2>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         1,
                                                         3>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       3>(i_1_1);
  Dune::XT::Functions::bind_combined_LocalizableFunction<YaspGrid2dType,
                                                         Dune::XT::Functions::internal::Combination::product,
                                                         1,
                                                         1,
                                                         1,
                                                         4>(m, yasp_grid_2d_id);
  Dune::XT::Functions::
      addbind_LocalizableFunctionInterface_combined_op<YaspGrid2dType,
                                                       Dune::XT::Functions::internal::Combination::product,
                                                       1,
                                                       1,
                                                       1,
                                                       4>(i_1_1);

  Dune::XT::Functions::bind_ConstantFunction<YaspGrid2dType, 1, 1>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ConstantFunction<YaspGrid2dType, 1, 2>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ConstantFunction<YaspGrid2dType, 1, 3>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ConstantFunction<YaspGrid2dType, 1, 4>(m, yasp_grid_2d_id);

  Dune::XT::Functions::bind_CheckerboardFunction<YaspGrid2dType, 1, 1>(m, yasp_grid_2d_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<YaspGrid2dType, 1, 2>(m, yasp_grid_2d_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<YaspGrid2dType, 1, 3>(m, yasp_grid_2d_id);
  //  Dune::XT::Functions::bind_CheckerboardFunction<YaspGrid2dType, 1, 4>(m, yasp_grid_2d_id);

  Dune::XT::Functions::bind_ExpressionFunction<YaspGrid2dType, 1, 1>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ExpressionFunction<YaspGrid2dType, 1, 2>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ExpressionFunction<YaspGrid2dType, 1, 3>(m, yasp_grid_2d_id);
  Dune::XT::Functions::bind_ExpressionFunction<YaspGrid2dType, 1, 4>(m, yasp_grid_2d_id);

  Dune::XT::Functions::bind_Spe10Model1Function<YaspGrid2dType, 1, 1>(m, yasp_grid_2d_id);

  return m.ptr();
}


#endif // HAVE_DUNE_PYBINDXI
