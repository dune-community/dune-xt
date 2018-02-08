// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/filesystem.hh>
#include <dune/xt/grid/grids.hh>

using namespace Dune;

template <class G, bool anything = true>
struct ExpectedResults;


struct Expected3dCubeResults
{
  static std::vector<size_t> local_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#else
    return std::vector<size_t>(27, 27);
#endif
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {{0, 1},  {1, 1},  {2, 1},  {3, 1},  {4, 1},  {5, 1},  {6, 1},  {7, 1},  {8, 1},
            {9, 1},  {10, 1}, {11, 1}, {12, 1}, {14, 1}, {15, 1}, {16, 1}, {17, 1}, {18, 1},
            {19, 1}, {20, 1}, {21, 1}, {22, 1}, {23, 1}, {24, 1}, {25, 1}, {26, 1}};
#else
    return {{0, 19},  {1, 15},  {2, 19},  {3, 15}, {4, 9},   {5, 15},  {6, 19},  {7, 15},  {8, 19},
            {9, 15},  {10, 9},  {11, 15}, {12, 9}, {14, 9},  {15, 15}, {16, 9},  {17, 15}, {18, 19},
            {19, 15}, {20, 19}, {21, 15}, {22, 9}, {23, 15}, {24, 19}, {25, 15}, {26, 19}};
#endif
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {{{1, 1}, {3, 1}, {9, 1}},
            {{0, 1}, {2, 1}, {4, 1}, {10, 1}},
            {{1, 1}, {5, 1}, {11, 1}},
            {{0, 1}, {4, 1}, {6, 1}, {12, 1}},
            {{1, 1}, {3, 1}, {5, 1}, {7, 1}, {13, 1}},
            {{2, 1}, {4, 1}, {8, 1}, {14, 1}},
            {{3, 1}, {7, 1}, {15, 1}},
            {{4, 1}, {6, 1}, {8, 1}, {16, 1}},
            {{5, 1}, {7, 1}, {17, 1}},
            {{0, 1}, {10, 1}, {12, 1}, {18, 1}},
            {{1, 1}, {9, 1}, {11, 1}, {13, 1}, {19, 1}},
            {{2, 1}, {10, 1}, {14, 1}, {20, 1}},
            {{3, 1}, {9, 1}, {13, 1}, {15, 1}, {21, 1}},
            {{4, 1}, {10, 1}, {12, 1}, {14, 1}, {16, 1}, {22, 1}},
            {{5, 1}, {11, 1}, {13, 1}, {17, 1}, {23, 1}},
            {{6, 1}, {12, 1}, {16, 1}, {24, 1}},
            {{7, 1}, {13, 1}, {15, 1}, {17, 1}, {25, 1}},
            {{8, 1}, {14, 1}, {16, 1}, {26, 1}},
            {{9, 1}, {19, 1}, {21, 1}},
            {{10, 1}, {18, 1}, {20, 1}, {22, 1}},
            {{11, 1}, {19, 1}, {23, 1}},
            {{12, 1}, {18, 1}, {22, 1}, {24, 1}},
            {{13, 1}, {19, 1}, {21, 1}, {23, 1}, {25, 1}},
            {{14, 1}, {20, 1}, {22, 1}, {26, 1}},
            {{15, 1}, {21, 1}, {25, 1}},
            {{16, 1}, {22, 1}, {24, 1}, {26, 1}},
            {{17, 1}, {23, 1}, {25, 1}}};
#else
    return {{{1, 9}, {3, 9}, {9, 9}},
            {{0, 9}, {2, 9}, {4, 9}, {10, 9}},
            {{1, 9}, {5, 9}, {11, 9}},
            {{0, 9}, {4, 9}, {6, 9}, {12, 9}},
            {{1, 9}, {3, 9}, {5, 9}, {7, 9}, {13, 9}},
            {{2, 9}, {4, 9}, {8, 9}, {14, 9}},
            {{3, 9}, {7, 9}, {15, 9}},
            {{4, 9}, {6, 9}, {8, 9}, {16, 9}},
            {{5, 9}, {7, 9}, {17, 9}},
            {{0, 9}, {10, 9}, {12, 9}, {18, 9}},
            {{1, 9}, {9, 9}, {11, 9}, {13, 9}, {19, 9}},
            {{2, 9}, {10, 9}, {14, 9}, {20, 9}},
            {{3, 9}, {9, 9}, {13, 9}, {15, 9}, {21, 9}},
            {{4, 9}, {10, 9}, {12, 9}, {14, 9}, {16, 9}, {22, 9}},
            {{5, 9}, {11, 9}, {13, 9}, {17, 9}, {23, 9}},
            {{6, 9}, {12, 9}, {16, 9}, {24, 9}},
            {{7, 9}, {13, 9}, {15, 9}, {17, 9}, {25, 9}},
            {{8, 9}, {14, 9}, {16, 9}, {26, 9}},
            {{9, 9}, {19, 9}, {21, 9}},
            {{10, 9}, {18, 9}, {20, 9}, {22, 9}},
            {{11, 9}, {19, 9}, {23, 9}},
            {{12, 9}, {18, 9}, {22, 9}, {24, 9}},
            {{13, 9}, {19, 9}, {21, 9}, {23, 9}, {25, 9}},
            {{14, 9}, {20, 9}, {22, 9}, {26, 9}},
            {{15, 9}, {21, 9}, {25, 9}},
            {{16, 9}, {22, 9}, {24, 9}, {26, 9}},
            {{17, 9}, {23, 9}, {25, 9}}};
#endif
  }
}; // Expected3dCubeResults


struct Expected3dSimplexResults
{
  static std::vector<size_t> local_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6};
#else
    return std::vector<size_t>(27, 162);
#endif
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {{0, 6},  {1, 4},  {2, 4},  {3, 4},  {4, 2},  {5, 3},  {6, 4},  {7, 3},  {8, 4},
            {9, 4},  {10, 2}, {11, 3}, {12, 2}, {14, 2}, {15, 3}, {16, 2}, {17, 4}, {18, 4},
            {19, 3}, {20, 4}, {21, 3}, {22, 2}, {23, 4}, {24, 4}, {25, 4}, {26, 6}};
#else
    return {{0, 54},  {1, 36},  {2, 48},  {3, 36},  {4, 18},  {5, 33},  {6, 48},  {7, 33},  {8, 48},
            {9, 36},  {10, 18}, {11, 33}, {12, 18}, {14, 18}, {15, 33}, {16, 18}, {17, 36}, {18, 48},
            {19, 33}, {20, 48}, {21, 33}, {22, 18}, {23, 36}, {24, 48}, {25, 36}, {26, 54}};
#endif
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
#if DXT_DISABLE_LARGE_TESTS
    return {{{1, 2}, {3, 2}, {9, 2}},
            {{0, 2}, {2, 2}, {4, 2}, {10, 2}},
            {{1, 2}, {5, 2}, {11, 2}},
            {{0, 2}, {4, 2}, {6, 2}, {12, 2}},
            {{1, 2}, {3, 2}, {5, 2}, {7, 2}, {13, 2}},
            {{2, 2}, {4, 2}, {8, 2}, {14, 2}},
            {{3, 2}, {7, 2}, {15, 2}},
            {{4, 2}, {6, 2}, {8, 2}, {16, 2}},
            {{5, 2}, {7, 2}, {17, 2}},
            {{0, 2}, {10, 2}, {12, 2}, {18, 2}},
            {{1, 2}, {9, 2}, {11, 2}, {13, 2}, {19, 2}},
            {{2, 2}, {10, 2}, {14, 2}, {20, 2}},
            {{3, 2}, {9, 2}, {13, 2}, {15, 2}, {21, 2}},
            {{4, 2}, {10, 2}, {12, 2}, {14, 2}, {16, 2}, {22, 2}},
            {{5, 2}, {11, 2}, {13, 2}, {17, 2}, {23, 2}},
            {{6, 2}, {12, 2}, {16, 2}, {24, 2}},
            {{7, 2}, {13, 2}, {15, 2}, {17, 2}, {25, 2}},
            {{8, 2}, {14, 2}, {16, 2}, {26, 2}},
            {{9, 2}, {19, 2}, {21, 2}},
            {{10, 2}, {18, 2}, {20, 2}, {22, 2}},
            {{11, 2}, {19, 2}, {23, 2}},
            {{12, 2}, {18, 2}, {22, 2}, {24, 2}},
            {{13, 2}, {19, 2}, {21, 2}, {23, 2}, {25, 2}},
            {{14, 2}, {20, 2}, {22, 2}, {26, 2}},
            {{15, 2}, {21, 2}, {25, 2}},
            {{16, 2}, {22, 2}, {24, 2}, {26, 2}},
            {{17, 2}, {23, 2}, {25, 2}}};
#else
    return {{{1, 18}, {3, 18}, {9, 18}},
            {{0, 18}, {2, 18}, {4, 18}, {10, 18}},
            {{1, 18}, {5, 18}, {11, 18}},
            {{0, 18}, {4, 18}, {6, 18}, {12, 18}},
            {{1, 18}, {3, 18}, {5, 18}, {7, 18}, {13, 18}},
            {{2, 18}, {4, 18}, {8, 18}, {14, 18}},
            {{3, 18}, {7, 18}, {15, 18}},
            {{4, 18}, {6, 18}, {8, 18}, {16, 18}},
            {{5, 18}, {7, 18}, {17, 18}},
            {{0, 18}, {10, 18}, {12, 18}, {18, 18}},
            {{1, 18}, {9, 18}, {11, 18}, {13, 18}, {19, 18}},
            {{2, 18}, {10, 18}, {14, 18}, {20, 18}},
            {{3, 18}, {9, 18}, {13, 18}, {15, 18}, {21, 18}},
            {{4, 18}, {10, 18}, {12, 18}, {14, 18}, {16, 18}, {22, 18}},
            {{5, 18}, {11, 18}, {13, 18}, {17, 18}, {23, 18}},
            {{6, 18}, {12, 18}, {16, 18}, {24, 18}},
            {{7, 18}, {13, 18}, {15, 18}, {17, 18}, {25, 18}},
            {{8, 18}, {14, 18}, {16, 18}, {26, 18}},
            {{9, 18}, {19, 18}, {21, 18}},
            {{10, 18}, {18, 18}, {20, 18}, {22, 18}},
            {{11, 18}, {19, 18}, {23, 18}},
            {{12, 18}, {18, 18}, {22, 18}, {24, 18}},
            {{13, 18}, {19, 18}, {21, 18}, {23, 18}, {25, 18}},
            {{14, 18}, {20, 18}, {22, 18}, {26, 18}},
            {{15, 18}, {21, 18}, {25, 18}},
            {{16, 18}, {22, 18}, {24, 18}, {26, 18}},
            {{17, 18}, {23, 18}, {25, 18}}};
#endif
  }
}; // Expected3dSimplexResults


template <bool anything>
struct ExpectedResults<YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>, anything> : public Expected3dCubeResults
{
  static std::string grid_name()
  {
    return "yasp_3d";
  }
};

#if HAVE_DUNE_ALUGRID

template <bool anything>
struct ExpectedResults<Dune::ALUGrid<3, 3, cube, nonconforming>, anything> : public Expected3dCubeResults
{
  static std::string grid_name()
  {
    return "alu_3d_cube";
  }
};

template <Dune::ALUGridRefinementType ref, bool anything>
struct ExpectedResults<Dune::ALUGrid<3, 3, simplex, ref>, anything> : public Expected3dSimplexResults
{
  static std::string grid_name()
  {
    return std::string("alu_3d_simplex_") + (ref == Dune::ALUGridRefinementType::conforming ? "" : "non")
           + "conforming";
  }
};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_DUNE_UGGRID

template <bool anything>
struct ExpectedResults<UGGrid<3>, anything> : public Expected3dSimplexResults
{
  static std::string grid_name()
  {
    return "ug_3d";
  }
};

#endif // HAVE_DUNE_UGGRID
#if HAVE_ALBERTA

template <bool anything>
struct ExpectedResults<AlbertaGrid<3, 3>, anything> : public Expected3dSimplexResults
{
  static std::string grid_name()
  {
    return "alberta_3d";
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
