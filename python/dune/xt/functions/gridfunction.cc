// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/grid-function.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace bindings {


template <class E, size_t r = 1, size_t rC = 1, class R = double>
class GridFunction
{
  using G = XT::Grid::extract_grid_t<E>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = Functions::GridFunction<E, r, rC, R>;
  using base_type = Functions::GridFunctionInterface<E, r, rC, R>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  template <size_t r_ = r, size_t rC_ = rC, bool anything = true>
  struct dim_dependent /*<r, rC, anything>, the general case, neither scalar nor vector nor square matrix*/
  {
    static void init(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const typename RangeTypeSelector<R, r, rC>::type&>(), "constant_matrix"_a);
      c.def(py::init<const FunctionInterface<d, r, rC, R>&>(), "matrix_function"_a, py::keep_alive<1, 2>());
      c.def(py::init<const GridFunctionInterface<E, r, rC, R>&>(), "matrix_grid_function"_a, py::keep_alive<1, 2>());
    } // ... init(...)

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // without dimRange
      m.def(FactoryName.c_str(),
            [](const GP&, const typename RangeTypeSelector<R, r, rC>::type& constant_matrix) {
              return type(constant_matrix);
            },
            "grid"_a,
            "constant_matrix"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, r, rC, R>& matrix_function) { return type(matrix_function); },
            "grid"_a,
            "matrix_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&, const GridFunctionInterface<E, r, rC, R>& matrix_grid_function) {
              return type(matrix_grid_function);
            },
            "grid"_a,
            "matrix_grid_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GridFunctionInterface<E, r, rC, R>& matrix_grid_function) { return type(matrix_grid_function); },
            "matrix_grid_function"_a,
            py::keep_alive<0, 1>());
      // and with dimRange
      m.def(FactoryName.c_str(),
            [](const GP&,
               const typename RangeTypeSelector<R, r, rC>::type& constant_matrix,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(constant_matrix);
            },
            "grid"_a,
            "constant_matrix"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FunctionInterface<d, r, rC, R>& matrix_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(matrix_function);
            },
            "grid"_a,
            "matrix_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, r, rC, R>& matrix_grid_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(matrix_grid_function);
            },
            "grid"_a,
            "matrix_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
    } // ... factory (...)
  }; // struct dim_dependent<r, rC, anything>, the general case, neither scalar nor square matrix

  template <bool anything>
  struct dim_dependent<1, 1, anything> /*the scalar case*/
  {
    static void init(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const R&>(), "constant_scalar"_a);
      c.def(py::init<const FieldVector<R, 1>&>(), "constant_vector_of_length_one"_a);
      c.def(py::init<const FunctionInterface<d, 1, 1, R>&>(), "scalar_function"_a, py::keep_alive<1, 2>());
      c.def(py::init<const GridFunctionInterface<E, 1, 1, R>&>(), "scalar_grid_function"_a, py::keep_alive<1, 2>());
    } // ... init(...)

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // without dimRange
      m.def(FactoryName.c_str(),
            [](const GP&, const R& constant_scalar) { return type(constant_scalar); },
            "grid"_a,
            "constant_scalar"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FieldVector<R, 1>& constant_vector_of_length_one) {
              return type(constant_vector_of_length_one);
            },
            "grid"_a,
            "constant_vector_of_length_one"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, 1, 1, R>& scalar_function) { return type(scalar_function); },
            "grid"_a,
            "scalar_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&, const GridFunctionInterface<E, 1, 1, R>& scalar_grid_function) {
              return type(scalar_grid_function);
            },
            "grid"_a,
            "scalar_grid_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GridFunctionInterface<E, 1, 1, R>& scalar_grid_function) { return type(scalar_grid_function); },
            "scalar_grid_function"_a,
            py::keep_alive<0, 1>());
      // and with dimRange, to distinguish from the square matrix case
      m.def(FactoryName.c_str(),
            [](const GP&, const R& constant_scalar, const Grid::bindings::Dimension<1>&) {
              return type(constant_scalar);
            },
            "grid"_a,
            "constant_scalar"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FieldVector<R, 1>& constant_vector_of_length_one, const Grid::bindings::Dimension<1>&) {
              return type(constant_vector_of_length_one);
            },
            "grid"_a,
            "constant_vector_of_length_one"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, 1, 1, R>& scalar_function, const Grid::bindings::Dimension<1>&) {
              return type(scalar_function);
            },
            "grid"_a,
            "scalar_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, 1, 1, R>& scalar_grid_function,
               const Grid::bindings::Dimension<1>&) { return type(scalar_grid_function); },
            "grid"_a,
            "scalar_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      // and with dimRange x dimRangeCols
      m.def(FactoryName.c_str(),
            [](const GP&,
               const R& constant_scalar,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(constant_scalar);
            },
            "grid"_a,
            "constant_scalar"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FieldVector<R, 1>& constant_vector_of_length_one,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(constant_vector_of_length_one);
            },
            "grid"_a,
            "constant_vector_of_length_one"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FunctionInterface<d, 1, 1, R>& scalar_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(scalar_function);
            },
            "grid"_a,
            "scalar_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, 1, 1, R>& scalar_grid_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(scalar_grid_function);
            },
            "grid"_a,
            "scalar_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
    } // ... factory (...)
  }; // struct dim_dependent<1, 1, anything> the scalar case

  template <size_t r_, bool anything>
  struct dim_dependent<r_, 1, anything> /*the non scalar vector case*/
  {
    static void init(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<const typename RangeTypeSelector<R, r, rC>::type&>(), "constant_vector"_a);
      c.def(py::init<const FunctionInterface<d, r, rC, R>&>(), "vector_function"_a, py::keep_alive<1, 2>());
      c.def(py::init<const GridFunctionInterface<E, r, rC, R>&>(), "vector_grid_function"_a, py::keep_alive<1, 2>());
    } // ... init(...)

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // without dimRange
      m.def(FactoryName.c_str(),
            [](const GP&, const typename RangeTypeSelector<R, r, rC>::type& constant_vector) {
              return type(constant_vector);
            },
            "grid"_a,
            "constant_vector"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, r, rC, R>& vector_function) { return type(vector_function); },
            "grid"_a,
            "vector_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&, const GridFunctionInterface<E, r, rC, R>& vector_grid_function) {
              return type(vector_grid_function);
            },
            "grid"_a,
            "vector_grid_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GridFunctionInterface<E, r, rC, R>& vector_grid_function) { return type(vector_grid_function); },
            "vector_grid_function"_a,
            py::keep_alive<0, 1>());
      // with dimRange
      m.def(FactoryName.c_str(),
            [](const GP&,
               const typename RangeTypeSelector<R, r, rC>::type& constant_vector,
               const Grid::bindings::Dimension<r>&) { return type(constant_vector); },
            "grid"_a,
            "constant_vector"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, r, rC, R>& vector_function, const Grid::bindings::Dimension<r>&) {
              return type(vector_function);
            },
            "grid"_a,
            "vector_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, r, rC, R>& vector_grid_function,
               const Grid::bindings::Dimension<r>&) { return type(vector_grid_function); },
            "grid"_a,
            "vector_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      // with dimRange x 1
      m.def(FactoryName.c_str(),
            [](const GP&,
               const typename RangeTypeSelector<R, r, rC>::type& constant_vector,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<1>>&) {
              return type(constant_vector);
            },
            "grid"_a,
            "constant_vector"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FunctionInterface<d, r, rC, R>& vector_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<1>>&) {
              return type(vector_function);
            },
            "grid"_a,
            "vector_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, r, rC, R>& vector_grid_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<1>>&) {
              return type(vector_grid_function);
            },
            "grid"_a,
            "vector_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
    } // ... factory (...)
  }; // struct dim_dependent<r, 1, anything>, the non scalar vector case

  template <size_t r_, bool anything>
  struct dim_dependent<r_, r_, anything> /*the non scalar square matrix case*/
  {
    static void init(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // we have the scalar variants ...
      dim_dependent<1, 1>::init(c);

      // .. and the matrix ones
      c.def(py::init<const FieldMatrix<R, r, rC>&>(), "constant_square_matrix"_a);
      c.def(py::init<const FunctionInterface<d, r, rC, R>&>(), "square_matrix_function"_a, py::keep_alive<1, 2>());
      c.def(py::init<const GridFunctionInterface<E, r, rC, R>&>(),
            "square_matrix_grid_function"_a,
            py::keep_alive<1, 2>());
    } // ... init(...)

    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      // we have the scalar variants ...
      dim_dependent<1, 1>::factory(m, FactoryName);

      // .. and the matrix ones, without dimRange
      m.def(FactoryName.c_str(),
            [](const GP&, const FieldMatrix<R, r, rC>& constant_square_matrix) { return type(constant_square_matrix); },
            "grid"_a,
            "constant_square_matrix"_a);
      m.def(FactoryName.c_str(),
            [](const GP&, const FunctionInterface<d, r, rC, R>& square_matrix_function) {
              return type(square_matrix_function);
            },
            "grid"_a,
            "square_matrix_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&, const GridFunctionInterface<E, r, rC, R>& square_matrix_grid_function) {
              return type(square_matrix_grid_function);
            },
            "grid"_a,
            "square_matrix_grid_function"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GridFunctionInterface<E, r, rC, R>& square_matrix_grid_function) {
              return type(square_matrix_grid_function);
            },
            "square_matrix_grid_function"_a,
            py::keep_alive<0, 1>());
      // .. and with dimRange
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FieldMatrix<R, r, rC>& constant_square_matrix,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(constant_square_matrix);
            },
            "grid"_a,
            "constant_square_matrix"_a,
            "dim_range"_a);
      m.def(FactoryName.c_str(),
            [](const GP&,
               const FunctionInterface<d, r, rC, R>& square_matrix_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(square_matrix_function);
            },
            "grid"_a,
            "square_matrix_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
      m.def(FactoryName.c_str(),
            [](const GP&,
               const GridFunctionInterface<E, r, rC, R>& square_matrix_grid_function,
               const std::pair<Grid::bindings::Dimension<r>, Grid::bindings::Dimension<rC>>&) {
              return type(square_matrix_grid_function);
            },
            "grid"_a,
            "square_matrix_grid_function"_a,
            "dim_range"_a,
            py::keep_alive<0, 2>());
    } // ... factory (...)
  }; // struct dim_dependent<r_, r_, anything> the square matrix (but not scalar) case

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "grid_function",
                         const std::string& grid_id = Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_to_" + Common::to_string(r);
    if (rC > 1)
      class_name += "x" + Common::to_string(rC);
    class_name += "d";
    if (!std::is_same<R, double>::value)
      class_name += "_" + Common::Typename<R>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), Common::to_camel_case(class_id).c_str());
    dim_dependent<>::init(c);

    dim_dependent<>::factory(m, Common::to_camel_case(class_id));
    return c;
  }
}; // class GridFunction


} // namespace bindings
} // namespace Functions
} // namespace XT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct GridFunction_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::XT::Functions::bindings::GridFunction<E, 1, 1>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 1, 2>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 1, 3>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 2, 1>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 2, 2>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 2, 3>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 3, 1>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 3, 2>::bind(m);
    Dune::XT::Functions::bindings::GridFunction<E, 3, 3>::bind(m);

    GridFunction_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct GridFunction_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_functions_gridfunction, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.functions._functions_gridfunction_interface_1d");
  py::module::import("dune.xt.functions._functions_gridfunction_interface_2d");
  py::module::import("dune.xt.functions._functions_gridfunction_interface_3d");

  GridFunction_for_all_grids<>::bind(m);
} // PYBIND11_MODULE(...)
