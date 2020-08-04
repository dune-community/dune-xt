// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_AS_GRID_FUNCTION_HH
#define PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_AS_GRID_FUNCTION_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/gridprovider.hh>

#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class G, size_t d, size_t r, size_t rC>
auto bind_FunctionAsGridFunctionWrapper(pybind11::module& m, const std::string& grid_id)
{
  static_assert(Grid::is_grid<G>::value);
  namespace py = pybind11;
  using namespace pybind11::literals;

  using E = typename G::template Codim<0>::Entity;
  using R = double;
  using I = GridFunctionInterface<E, r, rC, R>;
  using C = FunctionAsGridFunctionWrapper<E, r, rC, R>;

  const std::string classname = std::string("FunctionAsGridFunctionWrapper__" + grid_id + "_to_" + Common::to_string(r)
                                            + "x" + Common::to_string(rC));
  py::class_<C, I> c(m, classname.c_str(), classname.c_str());

  m.def(
      "function_to_grid_function",
      [](XT::Functions::FunctionInterface<d, r, rC, R>& func,
         const XT::Grid::GridProvider<G>& /*only_here_to_select_grid_type*/) { return std::make_unique<C>(func); },
      py::keep_alive<0, 1>(),
      "function"_a,
      "grid_provider"_a);

  return c;
} // ... bind_FunctionAsGridFunctionWrapper(...)


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_FUNCTIONS_FUNCTION_AS_GRID_FUNCTION_HH
