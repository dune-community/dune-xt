// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2021)
//   Tobias Leibner  (2020 - 2021)

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "dd_glued.hh"

#if HAVE_DUNE_GRID_GLUE

namespace Dune {
namespace XT {
namespace Grid {


template <>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                       YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>
{
  static int num_coarse_refinements()
  {
    return 0;
  }

  static int num_local_refinements()
  {
    return 2;
  }

  static std::string id()
  {
    return "2d_yaspgrid_yaspgrid";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
    return {4};
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, YaspGrid<2,
// EquidistantOffsetCoordinates<double, 2>>, anything>

#  if HAVE_DUNE_UGGRID || HAVE_UG

template <>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>>
{
  static int num_coarse_refinements()
  {
    return 0;
  }

  static int num_local_refinements()
  {
    return 2;
  }

  static std::string id()
  {
    return "2d_yaspgrid_uggrid";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
    return {4};
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>, anything>

#  endif // HAVE_DUNE_UGGRID || HAVE_UG
#  if HAVE_ALBERTA

template <>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, AlbertaGrid<2, 2>>
{
  static int num_coarse_refinements()
  {
    return 0;
  }

  static int num_local_refinements()
  {
    return 2;
  }

  static std::string id()
  {
    return "2d_yaspgrid_albertagrid";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
    return {2};
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, AlbertaGrid<2, 2>, anything>

#  endif // HAVE_ALBERTA


} // namespace Grid
} // namespace XT
} // namespace Dune


using namespace Dune;
using namespace Dune::XT::Grid;


// clang-format off
using GridTypes = ::testing::Types< std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                                     YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>
//#if HAVE_DUNE_ALUGRID // both the conforming as well as the nonconforming alugrid don't work
//                      , std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
//                                   ALUGrid<2, 2, simplex, conforming>>
//#endif
#if !HAVE_MPI && (HAVE_DUNE_UGGRID || HAVE_UG)
                        , std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>>
#endif
#if HAVE_ALBERTA
                        , std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, AlbertaGrid<2, 2>>
#endif
                        >; // clang-format on

TYPED_TEST_SUITE(GluedDdGridTest, GridTypes);
TYPED_TEST(GluedDdGridTest, setup_works)
{
  this->setup();
}
TYPED_TEST(GluedDdGridTest, visualize_is_callable)
{
  this->visualize_is_callable();
}
TYPED_TEST(GluedDdGridTest, couplings_are_of_correct_size)
{
  this->couplings_are_of_correct_size();
}
TYPED_TEST(GluedDdGridTest, local_grids_are_constructable)
{
  this->local_grids_are_constructable();
}
TYPED_TEST(GluedDdGridTest, walk_coupling_views)
{
  this->walk_coupling_views();
}

#endif // HAVE_DUNE_GRID_GLUE
