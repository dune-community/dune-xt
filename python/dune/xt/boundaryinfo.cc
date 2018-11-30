// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#include "config.h"

#include <string>
#include <vector>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/common/tuple.hh>

#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/boundaryinfo.bindings.hh>
#include <python/dune/xt/grid/available_types.hh>


using namespace Dune::XT::Common;
using namespace Dune::XT::Grid;
using InfoTypes = template_tuple<tplwrap<AllDirichletBoundaryInfo>,
                                 tplwrap<AllNeumannBoundaryInfo>,
                                 tplwrap<BoundarySegmentIndexBasedBoundaryInfo>,
                                 tplwrap<NormalBasedBoundaryInfo>>;
const std::array<std::string, 4> info_names{"all_dirichlet_boundary_info",
                                            "all_neumann_boundary_info",
                                            "boundary_segment_index_based_boundary_info",
                                            "normal_based_boundary_info"};
constexpr std::array<Layers, 6> layers{Layers::leaf,
                                       Layers::level,
                                       Layers::dd_subdomain,
                                       Layers::dd_subdomain_oversampled,
                                       Layers::dd_subdomain_boundary,
                                       Layers::dd_subdomain_coupling};

template <class Grid, Layers layer, class InfoTuple = InfoTypes>
struct bind_grid_layer_info
{
  static void bind(pybind11::module& m, size_t info_name_no = 0)
  {
    using I = extract_intersection_t<
        typename Dune::XT::Grid::
            Layer<Grid, layer, Dune::XT::Grid::Backends::view, Dune::XT::Grid::DD::SubdomainGrid<Grid>>::type>;
    using Info = typename InfoTuple::template head_type<I>;
    using binder = Dune::XT::Grid::bindings::BoundaryInfo<Info, Grid, layer>;
    binder::bind(m, info_names[info_name_no], layer_names[layer]);

    using Tail = typename InfoTuple::template tail_type<I>;
    bind_grid_layer_info<Grid, layer, Tail>::bind(m, ++info_name_no);
  }
};

template <class Grid, Layers l>
struct bind_grid_layer_info<Grid, l, null_template_tuple>
{
  static void bind(pybind11::module&, size_t)
  {
  }
};


template <class>
void bind_grid_layer(pybind11::module&, std::integral_constant<size_t, 0>)
{
}

template <class Grid, size_t counter>
void bind_grid_layer(pybind11::module& m, std::integral_constant<size_t, counter>)
{
  bind_grid_layer_info<Grid, layers[counter - 1]>::bind(m);
  bind_grid_layer<Grid>(m, std::integral_constant<size_t, counter - 1>());
}

template <class GridTuple = Dune::XT::Grid::bindings::AvailableTypes>
void bind_grid(pybind11::module& m)
{
  using Grid = typename GridTuple::head_type;
  bind_grid_layer<Grid>(m, std::integral_constant<size_t, layers.size()>());
  bind_grid<typename GridTuple::tail_type>(m);
}

template <>
void bind_grid<boost::tuples::null_type>(pybind11::module&)
{
}

PYBIND11_MODULE(_boundaryinfo, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  bind_grid(m);
}
