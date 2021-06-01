// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef PYTHON_DUNE_XT_GRID_TRAITS_HH
#define PYTHON_DUNE_XT_GRID_TRAITS_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/cast.h>

#include <dune/geometry/type.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/grid/dd/glued.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune::XT::Grid::bindings {


class Simplex
{};


class Cube
{};


class Pyramid
{};


class Prism
{};


template <size_t d>
class Dimension
{};


class LeafIntersection
{};


template <class MacroGridType, class MicroGridType>
class CouplingIntersection
{};


template <class MacroGridType, class MicroGridType>
class CouplingIntersectionBinder
{
public:
  using type = CouplingIntersection<MacroGridType, MicroGridType>;
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& macro_grid_id = grid_name<MacroGridType>::value(),
                         const std::string& micro_grid_id = grid_name<MicroGridType>::value(),
                         const std::string& class_id = "coupling_intersection")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + macro_grid_id + "_" + micro_grid_id;
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init(
        [](XT::Grid::DD::Glued<MacroGridType, MicroGridType, XT::Grid::Layers::leaf>&) { return new type(); }));

    // factory
    m.def(Common::to_camel_case(class_id).c_str(),
          [](XT::Grid::DD::Glued<MacroGridType, MicroGridType, XT::Grid::Layers::leaf>&) { return new type(); });

    return c;
  } // ... bind(...)
};


} // namespace Dune::XT::Grid::bindings


#endif // PYTHON_DUNE_XT_GRID_TRAITS_HH
