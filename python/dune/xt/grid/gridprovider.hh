// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_PBH
#define DUNE_XT_GRID_GRIDPROVIDER_PBH

#include <algorithm>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <dune/xt/common/parallel/mpi_comm_wrapper.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/entity.hh>

#include <dune/xt/grid/gridprovider.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class G>
void bind_GridProvider(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef GridProviderFactory<G> F;
  const std::string fac_id = std::string("GridProviderFactory__" + grid_id);
  py::class_<F> fac(m, fac_id.c_str(), fac_id.c_str());
  fac.def_static("available", &F::available);
  fac.def_static("default_config", &F::default_config);


  typedef GridProvider<G, Grid::none_t> C;

  py::class_<C> c(m, std::string("GridProvider__" + grid_id).c_str(), std::string("GridProvider__" + grid_id).c_str());

  c.def("max_level", &C::max_level);
  c.def("global_refine", &C::global_refine, "count"_a);
  c.def("visualize",
        [](const C& self, const std::string& filename, const Common::Configuration& boundary_info_cfg) {
          self.visualize(filename, boundary_info_cfg);
        },
        "filename"_a = C::static_id(),
        "boundary_info_cfg"_a = Common::Configuration());
  c.def_property_readonly("num_elements", [](const C& self) { return self.leaf_view().indexSet().size(0); });
  c.def_property_readonly("num_subdomains", [](const C& /*self*/) { return 1; });
  c.def("max_entity_diameter", [](const C& self) {
    double max_h = std::numeric_limits<double>::min();
    auto view = self.leaf_view();
    Walker<decltype(view)> walker(view);
    walker.append([&](const auto& entity) { max_h = std::max(max_h, entity_diameter(entity)); });
    walker.walk();
    return max_h;
  });
  c.def_property_readonly("grid_type", [grid_id](const C& /*self*/) { return grid_id; });
  c.def_property_readonly("dim", [](const C& /*self*/) { return C::dimDomain; });
} // ... bind_GridProvider(...)


template <class G>
void bind_DdSubdomainsGridProvider(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef GridProvider<G, DD::SubdomainGrid<G>> C;

  py::class_<C> c(m,
                  std::string("GridProviderWithDdSubdomains__" + grid_id).c_str(),
                  std::string("GridProviderWithDdSubdomains__" + grid_id).c_str());

  c.def("visualize",
        [](const C& self, const std::string& filename, const bool with_coupling) {
          self.visualize_dd(filename, with_coupling);
        },
        "filename"_a = C::static_id(),
        "with_coupling"_a = true);
  c.def_property_readonly("num_elements", [](const C& self) { return self.leaf_view().indexSet().size(0); });
  c.def("max_entity_diameter", [](const C& self) {
    double max_h = std::numeric_limits<double>::min();
    auto view = self.leaf_view();
    Walker<decltype(view)> walker(view);
    walker.append([&](const auto& entity) { max_h = std::max(max_h, entity_diameter(entity)); });
    walker.walk();
    return max_h;
  });
  c.def_property_readonly("num_subdomains", [](const C& self) { return self.dd_grid().size(); });
  c.def_property_readonly("grid_type", [grid_id](const C& /*self*/) { return grid_id; });
  c.def_property_readonly("dim", [](const C& /*self*/) { return C::dimDomain; });
  c.def("subdomain_on_boundary",
        [](const C& self, const ssize_t subdomain) {
          return self.dd_grid().boundary(XT::Common::numeric_cast<size_t>(subdomain));
        },
        "subdomain"_a);
  c.def("boundary_subdomains", [](const C& self) {
    std::vector<size_t> boundary_subdomains;
    for (size_t ss = 0; ss < self.dd_grid().size(); ++ss)
      if (self.dd_grid().boundary(ss))
        boundary_subdomains.emplace_back(ss);
    return boundary_subdomains;
  });
  c.def("neighboring_subdomains", [](const C& self, const ssize_t subdomain) {
    return self.dd_grid().neighborsOf(XT::Common::numeric_cast<size_t>(subdomain));
  });
  c.def("neighborhood_of",
        [](const C& self, const ssize_t subdomain) {
          const auto ss = XT::Common::numeric_cast<size_t>(subdomain);
          if (ss >= self.dd_grid().size())
            DUNE_THROW(XT::Common::Exceptions::index_out_of_range,
                       "subdomain: " << ss << "\n   self.num_subdomains: " << self.dd_grid().size());
          std::set<size_t> neighborhood;
          neighborhood.insert(ss);
          if (self.dd_grid().oversampling())
            for (auto&& entity : elements(self.dd_grid().local_grid_view(ss, /*oversampling=*/true)))
              neighborhood.insert(self.dd_grid().subdomainOf(entity));
          std::vector<size_t> result(neighborhood.begin(), neighborhood.end());
          std::sort(result.begin(), result.end());
          return result;
        },
        "subdomain"_a);
} // ... bind_DdSubdomainsGridProvider(...)


template <class G>
void bind_make_cube_grid(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.def(std::string("make_cube_grid__" + grid_id).c_str(),
        [](const Common::Configuration& cfg, Common::MPI_Comm_Wrapper mpi_comm) {
          return make_cube_grid<G>(cfg, mpi_comm.get());
        },
        "cfg"_a = cube_gridprovider_default_config(),
        "mpi_comm"_a = Common::MPI_Comm_Wrapper());

  m.def(std::string("make_cube_grid__" + grid_id).c_str(),
        [](const FieldVector<typename G::ctype, G::dimension>& lower_left,
           const FieldVector<typename G::ctype, G::dimension>& upper_right,
           const std::array<unsigned int, G::dimension>& num_elements,
           const unsigned int num_refinements,
           const std::array<unsigned int, G::dimension>& overlap_size,
           Common::MPI_Comm_Wrapper mpi_comm) {
          return make_cube_grid<G>(
              lower_left, upper_right, num_elements, num_refinements, overlap_size, mpi_comm.get());
        },
        "lower_left"_a,
        "upper_right"_a,
        "num_elements"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_dd_subdomains_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements")),
        "num_refinements"_a = cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
        "overlap_size"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_dd_subdomains_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size")),
        "mpi_comm"_a = Common::MPI_Comm_Wrapper());
} // ... bind_make_cube_grid(...)


template <class G>
void bind_make_cube_dd_subdomains_grid(pybind11::module& m, const std::string& grid_id)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.def(std::string("make_cube_dd_subdomains_grid__" + grid_id).c_str(),
        [](const FieldVector<typename G::ctype, G::dimension>& lower_left,
           const FieldVector<typename G::ctype, G::dimension>& upper_right,
           const std::array<unsigned int, G::dimension>& num_elements,
           const unsigned int num_refinements,
           const std::array<unsigned int, G::dimension>& overlap_size,
           const std::array<unsigned int, G::dimension> num_partitions,
           const size_t num_oversampling_layers,
           const size_t inner_boundary_segment_index) {
          return make_cube_dd_subdomains_grid<G>(lower_left,
                                                 upper_right,
                                                 num_elements,
                                                 num_refinements,
                                                 overlap_size,
                                                 num_partitions,
                                                 num_oversampling_layers,
                                                 inner_boundary_segment_index);
        },
        "lower_left"_a,
        "upper_right"_a,
        "num_elements"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_dd_subdomains_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements")),
        "num_refinements"_a = cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
        "overlap_size"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_dd_subdomains_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size")),
        "num_partitions"_a = XT::Common::make_array<unsigned int, G::dimension>(
            cube_dd_subdomains_gridprovider_default_config().template get<std::vector<unsigned int>>("num_partitions")),
        "num_oversampling_layers"_a =
            cube_dd_subdomains_gridprovider_default_config().template get<size_t>("num_refinements"),
        "inner_boundary_segment_index"_a =
            cube_dd_subdomains_gridprovider_default_config().template get<size_t>("inner_boundary_segment_index"));
} // ... bind_make_cube_dd_subdomains_grid(...)


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_PBH
