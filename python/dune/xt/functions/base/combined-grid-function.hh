// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTION_HH
#define PYTHON_DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTION_HH

#include <python/dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class G, class E, size_t r = 1, size_t rC = 1, class R = double>
class DifferenceGridFunction : public GridFunctionInterface<G, E, r, rC, R>
{
  using BaseType = GridFunctionInterface<G, E, r, rC, R>;

public:
  using base_type = typename BaseType::type;
  using type = Functions::DifferenceGridFunction<base_type, base_type>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "difference_grid_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = Common::to_camel_case(BaseType::class_name(grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());

    c.def(py::init<const base_type&, const base_type&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    BaseType::addbind_methods(c);

    return c;
  }
}; // class DifferenceGridFunction


template <class G, class E, size_t r = 1, size_t rC = 1, class R = double>
class SumGridFunction : public GridFunctionInterface<G, E, r, rC, R>
{
  using BaseType = GridFunctionInterface<G, E, r, rC, R>;

public:
  using base_type = typename BaseType::type;
  using type = Functions::SumGridFunction<base_type, base_type>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "sum_grid_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = Common::to_camel_case(BaseType::class_name(grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());

    c.def(py::init<const base_type&, const base_type&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    BaseType::addbind_methods(c);

    return c;
  }
}; // class SumGridFunction


template <class G, class E, class R = double>
class FractionGridFunction : public GridFunctionInterface<G, E, 1, 1, R>
{
  using BaseType = GridFunctionInterface<G, E, 1, 1, R>;

public:
  using base_type = typename BaseType::type;
  using type = Functions::FractionGridFunction<base_type, base_type>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "fraction_grid_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = Common::to_camel_case(BaseType::class_name(grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());

    c.def(py::init<const base_type&, const base_type&>(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>());

    BaseType::addbind_methods(c);

    return c;
  }
}; // class FractionGridFunction


template <class G, class E, size_t r = 1, size_t rC = 1, class R = double>
class ProductGridFunction : public GridFunctionInterface<G, E, r, rC, R>
{
  using BaseType = GridFunctionInterface<G, E, r, rC, R>;

public:
  using base_type = typename BaseType::type;
  using type = Functions::ProductGridFunction<Functions::GridFunctionInterface<E, 1, 1, R>, base_type>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "product_grid_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = Common::to_camel_case(BaseType::class_name(grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());

    c.def(py::init<const Functions::GridFunctionInterface<E, 1, 1, R>&, const base_type&>(),
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());

    BaseType::addbind_methods(c);

    return c;
  }
}; // class ProductGridFunction


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // PYTHON_DUNE_XT_FUNCTIONS_BASE_COMBINED_GRID_FUNCTION_HH
