// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2017 - 2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/filesystem.hh>
#include <dune/xt/grid/grids.hh>

using namespace Dune;

template <class G, bool anything = true>
struct ExpectedResults;

struct Expected1dResults
{
  static std::vector<size_t> local_sizes()
  {
    return {3, 3, 3};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 1}, {2, 1}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 1}}, {{0, 1}, {2, 1}}, {{1, 1}}};
  }
};


template <bool anything>
struct ExpectedResults<YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "yasp_1d";
  }
};

template <bool anything>
struct ExpectedResults<OneDGrid, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "oned_1d";
  }
};

#if HAVE_ALBERTA

template <bool anything>
struct ExpectedResults<AlbertaGrid<1, 1>, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "alberta_1d";
  }
};

#endif // HAVE_ALBERTA

#include "dd_subdomains_cube.hh"

TEST_F(CubeProviderTest, setup_works)
{
  this->setup();
}
TEST_F(CubeProviderTest, visualize_is_callable)
{
  this->visualize_is_callable(XT::Common::filename_only(::testing::internal::GetInjectableArgvs().at(0)));
}
TEST_F(CubeProviderTest, global_grid_part_entity_to_subdomain_mapping_is_correct)
{
  this->global_grid_part_entity_to_subdomain_mapping_is_correct();
}
TEST_F(CubeProviderTest, global_grid_view_entity_to_subdomain_mapping_is_correct)
{
  this->global_grid_view_entity_to_subdomain_mapping_is_correct();
}

TEST_F(CubeProviderTest, local_parts_are_of_correct_size)
{
  this->local_parts_are_of_correct_size();
}
TEST_F(CubeProviderTest, local_views_are_of_correct_size)
{
  this->local_views_are_of_correct_size();
}
TEST_F(CubeProviderTest, local_parts_are_indexed_consecutively)
{
  this->local_parts_are_indexed_consecutively();
}
TEST_F(CubeProviderTest, local_views_are_indexed_consecutively)
{
  this->local_views_are_indexed_consecutively();
}
TEST_F(CubeProviderTest, local_parts_report_correct_boundary_id)
{
  this->local_parts_report_correct_boundary_id();
}
TEST_F(CubeProviderTest, local_views_report_correct_boundary_id)
{
  this->local_views_report_correct_boundary_id();
}

TEST_F(CubeProviderTest, boundary_parts_are_of_correct_size)
{
  this->boundary_parts_are_of_correct_size();
}
TEST_F(CubeProviderTest, boundary_parts_are_indexed_consecutively)
{
  this->boundary_parts_are_indexed_consecutively();
}
TEST_F(CubeProviderTest, boundary_parts_contain_only_boundary_entities_and_intersections)
{
  this->boundary_parts_contain_only_boundary_entities_and_intersections();
}
TEST_F(CubeProviderTest, domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries)
{
  this->domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries();
}

TEST_F(CubeProviderTest, coupling_parts_are_of_correct_size)
{
  this->coupling_parts_are_of_correct_size();
}
TEST_F(CubeProviderTest, coupling_parts_are_indexed_consecutively)
{
  this->coupling_parts_are_indexed_consecutively();
}
TEST_F(CubeProviderTest, coupling_parts_contain_only_inner_entities_and_intersections)
{
  this->coupling_parts_contain_only_inner_entities_and_intersections();
}
TEST_F(CubeProviderTest, subdomain_connections_are_exactly_covered_by_couplings)
{
  this->subdomain_connections_are_exactly_covered_by_couplings();
}
