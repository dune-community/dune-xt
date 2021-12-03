// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <sstream>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/operators.h>
#include <dune/pybindxi/stl.h>
#include <dune/pybindxi/numpy.h>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/print.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class GV>
class Element
{
  using G = XT::Grid::extract_grid_t<GV>;
  using D = typename G::ctype;
  static constexpr size_t d = G::dimension;

public:
  using type = XT::Grid::extract_entity_t<GV>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "element",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11; // NOLINT(misc-unused-alias-decls)
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    // properties
    c.def_property_readonly("dimension", [](type& /*self*/) { return int(type::dimension); });
    c.def_property_readonly("codimension", [](type& /*self*/) { return int(type::codimension); });
    c.def_property_readonly("level", [](type& self) { return self.level(); });
    c.def_property_readonly("has_father", [](type& self) { return self.hasFather(); });
    c.def_property_readonly("is_leaf", [](type& self) { return self.isLeaf(); });
    c.def_property_readonly("is_regular", [](type& self) { return self.isRegular(); });
    c.def_property_readonly("is_new", [](type& self) { return self.isNew(); });
    c.def_property_readonly("might_vanish", [](type& self) { return self.mightVanish(); });
    c.def_property_readonly("has_boundary_intersections", [](type& self) { return self.hasBoundaryIntersections(); });
    c.def_property_readonly("affine", [](type& self) { return self.geometry().affine(); });
    c.def_property_readonly("volume", [](type& self) { return self.geometry().volume(); });
    c.def_property_readonly("center", [](type& self) {
      py::array_t<double> result(/*shape=*/d);
      auto access_to_result = result.mutable_unchecked<1>();
      const auto center = self.geometry().center();
      for (size_t ii = 0; ii < d; ++ii)
        access_to_result(ii) = center[ii];
      return result;
    });
    c.def_property_readonly("corners", [](type& self) {
      using XT::Common::numeric_cast;
      const auto num_corners = numeric_cast<size_t>(self.geometry().corners());
      py::array_t<double> result(/*shape=*/{num_corners, d});
      auto access_to_result = result.mutable_unchecked<2>();
      for (size_t cc = 0; cc < num_corners; ++cc) {
        const auto corner = self.geometry().corner(numeric_cast<int>(cc));
        for (size_t ii = 0; ii < d; ++ii)
          access_to_result(cc, ii) = corner[ii];
      }
      return result;
    });

    // special methods/operators
    c.def("__str__", [](type& self) {
      std::stringstream ss;
      ss << print(self);
      return ss.str();
    });
    c.def("__repr__", [](type& self) {
      std::stringstream ss;
      ss << repr(self);
      return ss.str();
    });
    c.def(
        "__eq__", [](type& self, const type& other) { return self == other; }, py::is_operator());
    c.def(
        "__neq__", [](type& self, const type& other) { return self != other; }, py::is_operator());

    // methods
    c.def(
        "sub_entities", [](type& self, const unsigned int codim) { return self.subEntities(codim); }, "codim"_a);
    c.def(
        "to_global",
        [](type& self, const py::array_t<double>& x_local) {
          py::array_t<double> result;
          if (x_local.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(x_local.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(0) = " << x_local.shape(0) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<1>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t dd = 0; dd < d; ++dd)
              x_local_dune[dd] = access_to_x_local(dd);
            const auto x_global = self.geometry().global(x_local_dune);
            result = py::array_t<double>(/*shape=*/d);
            auto access_to_result = result.mutable_unchecked<1>();
            for (size_t dd = 0; dd < d; ++dd)
              access_to_result(dd) = x_global[dd];
          } else if (x_local.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(x_local.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(1) = " << x_local.shape(1) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<2>();
            const size_t num_points = access_to_x_local.shape(0);
            result = py::array_t<double>(/*shape=*/{num_points, d});
            auto access_to_result = result.mutable_unchecked<2>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t ii = 0; ii < num_points; ++ii) {
              for (size_t dd = 0; dd < d; ++dd)
                x_local_dune[dd] = access_to_x_local(ii, dd);
              const auto x_global = self.geometry().global(x_local_dune);
              for (size_t dd = 0; dd < d; ++dd)
                access_to_result(ii, dd) = x_global[dd];
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "x_local.ndim() = " << x_local.ndim()
                                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "x_local"_a);
    c.def(
        "to_local",
        [](type& self, const py::array_t<double>& x_global) {
          py::array_t<double> result;
          if (x_global.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(x_global.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_global.shape(0) = " << x_global.shape(0) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_global = x_global.unchecked<1>();
            XT::Common::FieldVector<double, d> x_global_dune;
            for (size_t dd = 0; dd < d; ++dd)
              x_global_dune[dd] = access_to_x_global(dd);
            const auto x_local = self.geometry().local(x_global_dune);
            result = py::array_t<double>(/*shape=*/d);
            auto access_to_result = result.mutable_unchecked<1>();
            for (size_t dd = 0; dd < d; ++dd)
              access_to_result(dd) = x_local[dd];
          } else if (x_global.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(x_global.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_global.shape(1) = " << x_global.shape(1) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_global = x_global.unchecked<2>();
            const size_t num_points = access_to_x_global.shape(0);
            result = py::array_t<double>(/*shape=*/{num_points, d});
            auto access_to_result = result.mutable_unchecked<2>();
            XT::Common::FieldVector<double, d> x_global_dune;
            for (size_t ii = 0; ii < num_points; ++ii) {
              for (size_t dd = 0; dd < d; ++dd)
                x_global_dune[dd] = access_to_x_global(ii, dd);
              const auto x_local = self.geometry().local(x_global_dune);
              for (size_t dd = 0; dd < d; ++dd)
                access_to_result(ii, dd) = x_local[dd];
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "x_global.ndim() = " << x_global.ndim()
                                            << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "x_global"_a);
    c.def(
        "integration_element",
        [](type& self, const py::array_t<double>& x_local) {
          py::array_t<double> result;
          if (x_local.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(x_local.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(0) = " << x_local.shape(0) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<1>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t dd = 0; dd < d; ++dd)
              x_local_dune[dd] = access_to_x_local(dd);
            result = py::array_t<double>(/*shape=*/1);
            auto access_to_result = result.mutable_unchecked<1>();
            access_to_result(0) = self.geometry().integrationElement(x_local_dune);
          } else if (x_local.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(x_local.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(1) = " << x_local.shape(1) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<2>();
            const size_t num_points = access_to_x_local.shape(0);
            result = py::array_t<double>(/*shape=*/{num_points, size_t(1)});
            auto access_to_result = result.mutable_unchecked<2>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t ii = 0; ii < num_points; ++ii) {
              for (size_t dd = 0; dd < d; ++dd)
                x_local_dune[dd] = access_to_x_local(ii, dd);
              access_to_result(ii, 0) = self.geometry().integrationElement(x_local_dune);
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "x_local.ndim() = " << x_local.ndim()
                                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "x_local"_a);
    c.def(
        "jacobian_transposed",
        [](type& self, const py::array_t<double>& x_local) {
          py::array_t<double> result;
          if (x_local.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(x_local.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(0) = " << x_local.shape(0) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<1>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t dd = 0; dd < d; ++dd)
              x_local_dune[dd] = access_to_x_local(dd);
            const auto jac = self.geometry().jacobianTransposed(x_local_dune);
            result = py::array_t<double>(/*shape=*/{d, d});
            auto access_to_result = result.mutable_unchecked<2>();
            for (size_t ii = 0; ii < d; ++ii)
              for (size_t jj = 0; jj < d; ++jj)
                access_to_result(ii, jj) = jac[ii][jj];
          } else if (x_local.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(x_local.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(1) = " << x_local.shape(1) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<2>();
            const size_t num_points = access_to_x_local.shape(0);
            result = py::array_t<double>(/*shape=*/{num_points, d, d});
            auto access_to_result = result.mutable_unchecked<3>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t pp = 0; pp < num_points; ++pp) {
              for (size_t dd = 0; dd < d; ++dd)
                x_local_dune[dd] = access_to_x_local(pp, dd);
              const auto jac = self.geometry().jacobianTransposed(x_local_dune);
              for (size_t ii = 0; ii < d; ++ii)
                for (size_t jj = 0; jj < d; ++jj)
                  access_to_result(pp, ii, jj) = jac[ii][jj];
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "x_local.ndim() = " << x_local.ndim()
                                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "x_local"_a);
    c.def(
        "jacobian_inverse_transposed",
        [](type& self, const py::array_t<double>& x_local) {
          py::array_t<double> result;
          if (x_local.ndim() == 1) {
            // a single point
            DUNE_THROW_IF(x_local.shape(0) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(0) = " << x_local.shape(0) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<1>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t dd = 0; dd < d; ++dd)
              x_local_dune[dd] = access_to_x_local(dd);
            const auto jac_inv = self.geometry().jacobianInverseTransposed(x_local_dune);
            result = py::array_t<double>(/*shape=*/{d, d});
            auto access_to_result = result.mutable_unchecked<2>();
            for (size_t ii = 0; ii < d; ++ii)
              for (size_t jj = 0; jj < d; ++jj)
                access_to_result(ii, jj) = jac_inv[ii][jj];
          } else if (x_local.ndim() == 2) {
            // a list of points
            DUNE_THROW_IF(x_local.shape(1) != d,
                          XT::Common::Exceptions::shapes_do_not_match,
                          "x_local.shape(1) = " << x_local.shape(1) << " (has to be " << size_t(d) << ")!");
            const auto& access_to_x_local = x_local.unchecked<2>();
            const size_t num_points = access_to_x_local.shape(0);
            result = py::array_t<double>(/*shape=*/{num_points, d, d});
            auto access_to_result = result.mutable_unchecked<3>();
            XT::Common::FieldVector<double, d> x_local_dune;
            for (size_t pp = 0; pp < num_points; ++pp) {
              for (size_t dd = 0; dd < d; ++dd)
                x_local_dune[dd] = access_to_x_local(pp, dd);
              const auto jac_inv = self.geometry().jacobianInverseTransposed(x_local_dune);
              for (size_t ii = 0; ii < d; ++ii)
                for (size_t jj = 0; jj < d; ++jj)
                  access_to_result(pp, ii, jj) = jac_inv[ii][jj];
            }
          } else
            DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                       "x_local.ndim() = " << x_local.ndim()
                                           << " (has to be 1 for a single point, 2 for a list of points)!");
          return result;
        },
        "x_local"_a);

    return c;
  } // ... bind(...)
}; // namespace bindings


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct Element_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    Dune::XT::Grid::bindings::Element<GV>::bind(m);
    Element_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct Element_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_grid_element, m)
{
  namespace py = pybind11;
  using namespace Dune;

  py::module::import("dune.xt.common");

  Element_for_all_grids<>::bind(m);
}
