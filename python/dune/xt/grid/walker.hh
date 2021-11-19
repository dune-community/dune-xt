// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2021)

#ifndef PYTHON_DUNE_XT_GRID_WALKER_HH
#define PYTHON_DUNE_XT_GRID_WALKER_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/functional.h>

#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/walker.hh>
#include <xt/common/python.hh>


namespace Dune::XT::Grid::bindings {


template <class GV>
class Walker
{
  using G = typename GV::Grid;
  static_assert(is_grid<G>::value);

public:
  using type = Grid::Walker<GV>;
  using base_type = Grid::ElementAndIntersectionFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  using E = typename type::ElementType;
  using I = typename type::IntersectionType;

public:
  template <class T, typename... options>
  static void addbind_methods(pybind11::class_<T, options...>& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "append",
        [](T& self, ElementFunctor<GV>& functor, const ElementFilter<GV>& filter) { self.append(functor, filter); },
        "element_functor"_a,
        "element_filter"_a = ApplyOn::AllElements<GV>());
    c.def(
        "append",
        [](T& self, std::function<void(const E&)> generic_element_function, const ElementFilter<GV>& filter) {
          self.append(/*prepare=*/[]() {}, /*apply_local=*/generic_element_function, /*finalize=*/[]() {}, filter);
        },
        "generic_element_function"_a,
        "element_filter"_a = ApplyOn::AllElements<GV>());
    c.def(
        "append",
        [](T& self, IntersectionFunctor<GV>& functor, const IntersectionFilter<GV>& filter) {
          self.append(functor, filter);
        },
        "intersection_functor"_a,
        "intersection_filter"_a = ApplyOn::AllIntersections<GV>());
    c.def(
        "append",
        [](T& self,
           std::function<void(const I&, const E&, const E&)> generic_intersection_function,
           const IntersectionFilter<GV>& filter) {
          self.append(/*prepare=*/[]() {}, /*apply_local=*/generic_intersection_function, /*finalize=*/[]() {}, filter);
        },
        "generic_intersection_function"_a,
        "intersection_filter"_a = ApplyOn::AllIntersections<GV>());
    c.def(
        "append",
        [](T& self,
           ElementAndIntersectionFunctor<GV>& functor,
           const IntersectionFilter<GV>& intersection_filter,
           const ElementFilter<GV>& element_filter) { self.append(functor, intersection_filter, element_filter); },
        "element_and_intersection_functor"_a,
        "intersection_filter"_a = ApplyOn::AllIntersections<GV>(),
        "element_filter"_a = ApplyOn::AllElements<GV>());
    c.def(
        "walk",
        [](T& self, const bool thread_parallel = false, const bool clear_functors_after_walk = true) {
          if (!HAVE_TBB && thread_parallel) {
            Common::bindings::warning(PyExc_RuntimeWarning, "Requested thread parallel walk without TBB support");
          }
          self.walk(thread_parallel, clear_functors_after_walk);
        },
        "thread_parallel"_a = false,
        "clear_functors_after_walk"_a = true,
        py::call_guard<py::gil_scoped_release>());
  } // ... addbind_methods(...)

  static bound_type bind_leaf(pybind11::module& m,
                              const std::string& grid_id = grid_name<G>::value(),
                              const std::string& layer_id = "",
                              const std::string& class_id = "Walker")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def(py::init([](GridProvider<G>& grid_provider) { return new type(grid_provider.leaf_view()); }),
          "grid_provider"_a);

    addbind_methods(c);

    return c;
  } // ... bind_leaf(...)

  static bound_type bind_coupling(pybind11::module& m,
                                  const std::string& grid_id = grid_name<G>::value(),
                                  const std::string& layer_id = "",
                                  const std::string& class_id = "Walker")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    if (!layer_id.empty())
      ClassName += "_" + layer_id;
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def(py::init([](CouplingGridProvider<GV>& grid_provider) { return new type(grid_provider.coupling_view()); }),
          "grid_provider"_a);

    addbind_methods(c);

    return c;
  } // ... bind_coupling(...)


  static void bind_leaf_factory(pybind11::module& m, const std::string& class_id = "Walker")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](GridProvider<G>& grid_provider) { return new type(grid_provider.leaf_view()); },
        "grid_provider"_a,
        py::keep_alive<0, 1>());
  } // ... bind_leaf_factory(...)

  static void bind_coupling_factory(pybind11::module& m, const std::string& class_id = "Walker")
  {
    using namespace pybind11::literals;
    m.def(
        Common::to_camel_case(class_id).c_str(),
        [](CouplingGridProvider<GV>& coupling_grid_provider) {
          return new type(coupling_grid_provider.coupling_view());
        },
        "coupling_grid_provider"_a);
  } // ... bind_coupling_factory(...)
}; // class Walker


} // namespace Dune::XT::Grid::bindings

#endif // PYTHON_DUNE_XT_GRID_WALKER_HH
