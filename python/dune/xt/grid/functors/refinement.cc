// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/xt/grid/functors/refinement.hh>

#include "interfaces.hh"


namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
class MaximumEntityVolumeRefineFunctorFunctor
{
  static_assert(is_grid<G>::value, "");
  using GV = typename G::LeafGridView;
  using I = extract_intersection_t<GV>;

public:
  using type = Grid::MaximumEntityVolumeRefineFunctorFunctor<GV>;
  using base_type = Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "maximum_element_volume_refine_functor",
                         const std::string& grid_id = grid_name<G>::value())
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    auto ClassId = Common::to_camel_case(class_id);
    auto ClassName = Common::to_camel_case(class_id + "_" + grid_id);
    bound_type c(m, ClassName.c_str(), std::string(ClassId + "( " + grid_id + " variant)").c_str());
    c.def(py::init([](GridProvider<G>& grid_provider, const double& volume) {
            return std::make_unique<type>(grid_provider.grid(), volume, 1.);
          }),
          "grid_provider"_a,
          "volume"_a,
          py::keep_alive<0, 1>());
    c.def("__repr__", [ClassId](type&) { return ClassId + "(grid_provider=\?\?\?, volume=\?\?\?)"; });
    c.def_property_readonly("result", [](const type& self) { return self.result(); });

    m.def(
        ClassId.c_str(),
        [](GridProvider<G>& grid_provider, const double& volume) {
          return std::make_unique<type>(grid_provider.grid(), volume, 1.);
        },
        "grid_provider"_a,
        "volume"_a,
        py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)
}; // class MaximumEntityVolumeRefineFunctorFunctor


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MaximumEntityVolumeRefineFunctorFunctor_for_all_grids
{
  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::MaximumEntityVolumeRefineFunctorFunctor<typename GridTypes::head_type>::bind(m);
    MaximumEntityVolumeRefineFunctorFunctor_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct MaximumEntityVolumeRefineFunctorFunctor_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_functors_refinement, m)
{
  namespace py = pybind11;

  py::module::import("dune.xt.grid._grid_boundaryinfo_interfaces");
  py::module::import("dune.xt.grid._grid_boundaryinfo_types");
  py::module::import("dune.xt.grid._grid_gridprovider_provider");
  py::module::import("dune.xt.grid._grid_functors_interfaces");

  MaximumEntityVolumeRefineFunctorFunctor_for_all_grids<>::bind(m);
}
