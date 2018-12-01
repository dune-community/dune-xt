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


struct Expected2dSimplexNonconformingResults
{
  static std::vector<size_t> local_sizes()
  {
    return {18, 18, 18, 18, 18, 18, 18, 18, 18};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 6}, {1, 3}, {2, 5}, {3, 3}, {5, 3}, {6, 5}, {7, 3}, {8, 6}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 3}, {3, 3}},
            {{0, 3}, {2, 3}, {4, 3}},
            {{1, 3}, {5, 3}},
            {{0, 3}, {4, 3}, {6, 3}},
            {{1, 3}, {3, 3}, {5, 3}, {7, 3}},
            {{2, 3}, {4, 3}, {8, 3}},
            {{3, 3}, {7, 3}},
            {{4, 3}, {6, 3}, {8, 3}},
            {{5, 3}, {7, 3}}};
  }
}; // Expected2dSimplexNonconformingResults


struct Expected2dCubeResults
{
  static std::vector<size_t> local_sizes()
  {
    return {9, 9, 9, 9, 9, 9, 9, 9, 9};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 5}, {1, 3}, {2, 5}, {3, 3}, {5, 3}, {6, 5}, {7, 3}, {8, 5}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 3}, {3, 3}},
            {{0, 3}, {2, 3}, {4, 3}},
            {{1, 3}, {5, 3}},
            {{0, 3}, {4, 3}, {6, 3}},
            {{1, 3}, {3, 3}, {5, 3}, {7, 3}},
            {{2, 3}, {4, 3}, {8, 3}},
            {{3, 3}, {7, 3}},
            {{4, 3}, {6, 3}, {8, 3}},
            {{5, 3}, {7, 3}}};
  }
}; // Expected2dCubeResults


template <bool anything>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, anything> : public Expected2dCubeResults
{
  static std::string grid_name()
  {
    return "yasp_2d";
  }
};

#if HAVE_DUNE_ALUGRID

template <bool anything>
struct ExpectedResults<Dune::ALUGrid<2, 2, cube, nonconforming>, anything> : public Expected2dCubeResults
{
  static std::string grid_name()
  {
    return "alu_2d_cube";
  }
};

template <bool anything>
struct ExpectedResults<Dune::ALUGrid<2, 2, simplex, conforming>, anything>
{
  static std::string grid_name()
  {
    return "alu_2d_simplex_conforming";
  }

  static std::vector<size_t> local_sizes()
  {
    return {36, 36, 36, 36, 36, 36, 36, 36, 36};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 6}, {1, 3}, {2, 6}, {3, 3}, {5, 3}, {6, 6}, {7, 3}, {8, 6}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 3}, {3, 3}},
            {{0, 3}, {2, 3}, {4, 3}},
            {{1, 3}, {5, 3}},
            {{0, 3}, {4, 3}, {6, 3}},
            {{1, 3}, {3, 3}, {5, 3}, {7, 3}},
            {{2, 3}, {4, 3}, {8, 3}},
            {{3, 3}, {7, 3}},
            {{4, 3}, {6, 3}, {8, 3}},
            {{5, 3}, {7, 3}}};
  }
}; // ExpectedResults<Dune::ALUGrid<2, 2, simplex, conforming>, ...>

template <bool anything>
struct ExpectedResults<Dune::ALUGrid<2, 2, simplex, nonconforming>, anything>
  : public Expected2dSimplexNonconformingResults
{
  static std::string grid_name()
  {
    return "alu_2d_simplex_nonconforming";
  }
};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_DUNE_UGGRID

template <bool anything>
struct ExpectedResults<UGGrid<2>, anything> : public Expected2dSimplexNonconformingResults
{
  static std::string grid_name()
  {
    return "ug_2d_simplex";
  }
};

#endif // HAVE_DUNE_UGGRID
#if HAVE_ALBERTA

template <bool anything>
struct ExpectedResults<AlbertaGrid<2, 2>, anything> : public Expected2dSimplexNonconformingResults
{
  static std::string grid_name()
  {
    return "alberta_2d";
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
