// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_GRID_TEST_DD_GLUED_HH
#define DUNE_XT_GRID_TEST_DD_GLUED_HH
#if HAVE_DUNE_GRID_GLUE

#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/dd/glued.hh>


template <class T>
std::string convert_to_initializer_list_str(const std::set<T>& results)
{
  std::stringstream out;
  if (results.size() == 0)
    out << "{}";
  else if (results.size() == 1)
    out << "{" << *results.begin() << "}";
  else {
    auto iterator = results.begin();
    out << "{" << *iterator;
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ", " << *iterator;
    }
    out << "}";
  }
  return out.str();
}

template <class L, class R>
std::string convert_to_initializer_list_str(const std::pair<L, R>& results)
{
  std::stringstream out;
  out << "{" << results.first << ", " << results.second << "}";
  return out.str();
}

template <class F, class S>
std::string convert_to_initializer_list_str(const std::map<F, S>& results)
{
  std::stringstream out;
  if (results.size() == 0)
    out << "{}" << std::endl;
  else if (results.size() == 1)
    out << "{{" << convert_to_initializer_list_str(results.begin()->first) << ", " << results.begin()->second << "}}";
  else {
    auto iterator = results.begin();
    out << "{{" << convert_to_initializer_list_str(iterator->first) << ", " << iterator->second << "}";
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ",\n {" << convert_to_initializer_list_str(iterator->first) << ", " << iterator->second << "}";
    }
    out << "}";
  }
  return out.str();
}


namespace Dune {
namespace XT {
namespace Grid {


template <class M, class L, bool anything = true>
struct ExpectedResults
{
  static_assert(AlwaysFalse<M>::value, "Please add me for this grid!");
};


/// \note assumes that all macro entities contain local grids of same refiment levels
template <class GridTuple>
struct GluedDdGridTest : public ::testing::Test
{
  typedef typename std::tuple_element<0, GridTuple>::type MacroGridType;
  typedef typename std::tuple_element<1, GridTuple>::type LocalGridType;
  typedef ExpectedResults<MacroGridType, LocalGridType> Expectations;

  template <class G, bool anything = true>
  struct get_local_layer
  {
    static const constexpr Layers type = Layers::level;
  };

#if HAVE_ALBERTA

  template <int d, int dw, bool anything>
  struct get_local_layer<AlbertaGrid<d, dw>, anything>
  {
    static const constexpr Layers type = Layers::leaf;
  };

#endif

  static const constexpr Layers local_layer = get_local_layer<LocalGridType>::type;

  void setup()
  {
    if (!macro_grid_)
      macro_grid_ = std::make_unique<GridProvider<MacroGridType>>(
          make_cube_grid<MacroGridType>(0., 1., 3, Expectations::num_coarse_refinements()));
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    if (!dd_grid_)
      dd_grid_ = std::make_unique<DD::Glued<MacroGridType, LocalGridType, local_layer>>(
          *macro_grid_,
          Expectations::num_local_refinements(),
          /*prepare_glues=*/false,
          /*allow_for_broken_orientation_of_coupling_intersections=*/true);
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    for (auto&& macro_entity : Dune::elements(dd_grid_->macro_grid_view())) {
      EXPECT_EQ(dd_grid_->max_local_level(macro_entity),
                (local_layer == Layers::level) ? Expectations::num_local_refinements() : -1);
    }
  } // ... setup()

  void couplings_are_of_correct_size()
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";

    const auto& macro_grid_view = dd_grid_->macro_grid_view();
    for (auto&& macro_entity : Dune::elements(macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection : Dune::intersections(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor = macro_intersection.outside();
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling =
              dd_grid_->coupling(macro_entity,
                                 (local_layer == Layers::level) ? 1 : -1,
                                 macro_neighbor,
                                 (local_layer == Layers::level) ? Expectations::num_local_refinements() : -1,
                                 /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          EXPECT_EQ(Expectations::num_local_couplings_intersections().count(coupling.size()), 1)
              << "entity: " << entity_index << "\n"
              << "neighbor: " << neighbor_index << "\n"
              << "expected num_local_couplings_intersections: "
              << convert_to_initializer_list_str(Expectations::num_local_couplings_intersections())
              << "\nactual num_local_couplings_intersections: " << coupling.size();
        }
      }
    }
  } // ... couplings_are_of_correct_size(...)

  void visualize_is_callable()
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";

    dd_grid_->visualize(Expectations::id());
  } // ... visualize_is_callable(...)

  void check_intersection_orientation_for_equal_levels(const bool expect_failure = Expectations::failure_for_equal())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    if (local_layer == Layers::level) {
      for (int level = 0; level <= dd_grid_->max_local_level(0); ++level)
        failure += check_intersections_for_levels(level, level, expect_failure);
    } else {
      failure += check_intersections_for_levels(-1, -1, expect_failure);
    }

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_equal_levels(...)

  void check_intersection_orientation_for_higher_neighbor_levels(
      const bool expect_failure = Expectations::failure_for_higher())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    if (local_layer == Layers::level) {
      for (int entity_level = 0; entity_level <= dd_grid_->max_local_level(0); ++entity_level)
        for (int neighbor_level = entity_level; neighbor_level <= dd_grid_->max_local_level(0); ++neighbor_level)
          failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);
    }

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_higher_neighbor_levels(...)

  void check_intersection_orientation_for_lower_or_equal_neighbor_levels(
      const bool expect_failure = Expectations::failure_for_lower_or_equal())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(dd_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    if (local_layer == Layers::level) {
      for (int entity_level = 0; entity_level <= dd_grid_->max_local_level(0); ++entity_level)
        for (int neighbor_level = 0; neighbor_level <= entity_level; ++neighbor_level)
          failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);
    }

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_lower_or_equal_neighbor_levels(...)

  size_t check_intersections_for_levels(const int entity_level, const int neighbor_level, const bool expect_failure)
  {
    size_t failure = 0;
    const auto actual_num_wrongly_oriented_intersections = count_wrong_intersections(entity_level, neighbor_level);
    if (expect_failure) {
      const auto expected_results = Expectations::results();
      const auto search_for_levels_in_expected_results =
          expected_results.find(std::make_pair(entity_level, neighbor_level));
      EXPECT_NE(search_for_levels_in_expected_results, expected_results.end())
          << "missing expected results for entity and neighbor level "
          << convert_to_initializer_list_str(std::make_pair(entity_level, neighbor_level)) << ".\n"
          << "-> PLEASE ADD A RECORD TO\n"
          << "   ExpectedResults<" << Common::Typename<MacroGridType>::value() << ", "
          << Common::Typename<LocalGridType>::value() << ">!";
      if (search_for_levels_in_expected_results == expected_results.end()) {
        DUNE_THROW(InvalidStateException,
                   "Cannot use ASSERT_EQ above, so we need to exit this way.\n\n"
                       << "The actual numbers of broken intersections are\n"
                       << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()));
      }
      const auto expected_num_wrongly_oriented_intersections = search_for_levels_in_expected_results->second;
      if (expected_num_wrongly_oriented_intersections != actual_num_wrongly_oriented_intersections)
        ++failure;
      EXPECT_EQ(expected_num_wrongly_oriented_intersections, actual_num_wrongly_oriented_intersections)
          << "\nTHIS IS A GOOD THING, AN ACTUAL NUMBER OF FAILURES WHICH IS LOWER THAN THE EXPECTED NUMBER OF "
          << "FAILURES IS AN IMPROVEMENT!\n"
          << "-> BE HAPPY AND UPDATE THE RECORD IN!\n"
          << "   ExpectedResults<" << Common::Typename<MacroGridType>::value() << ", "
          << Common::Typename<LocalGridType>::value() << ">!\n"
          << "IF THE UPDATED RECORDS DO NOT INDICATE FAILURES ANYMORE, UPDATE THE TESTS!\n";
    } else {
      if (actual_num_wrongly_oriented_intersections != 0)
        ++failure;
      EXPECT_EQ(actual_num_wrongly_oriented_intersections, 0);
    }
    return failure;
  } // ... check_intersections(...)

  std::map<std::pair<ssize_t, ssize_t>, size_t> count_wrong_intersections_on_all_levels()
  {
    setup();
    if (!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if (!dd_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- s.a.

    std::map<std::pair<ssize_t, ssize_t>, size_t> results;
    if (local_layer == Layers::level) {
      for (int entity_level = 0; entity_level <= dd_grid_->max_local_level(0); ++entity_level)
        for (int neighbor_level = 0; neighbor_level <= dd_grid_->max_local_level(0); ++neighbor_level)
          results[std::make_pair(entity_level, neighbor_level)] =
              count_wrong_intersections(entity_level, neighbor_level);
    } else {
      results[std::make_pair(-1, -1)] = count_wrong_intersections(-1, -1);
    }
    return results;
  }

  size_t count_wrong_intersections(const int entity_level, const int neighbor_level)
  {
    setup();
    if (!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if (!dd_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- s.a.

    const auto& macro_grid_view = dd_grid_->macro_grid_view();
    size_t failures = 0;
    for (auto&& macro_entity : Dune::elements(macro_grid_view)) {
      for (auto&& macro_intersection : Dune::intersections(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor = macro_intersection.outside();
          //          const auto local_grid_view = dd_grid_->local_grid(macro_entity).level_view(entity_level);
          const auto& coupling_glue =
              dd_grid_->coupling(macro_entity,
                                 entity_level,
                                 macro_neighbor,
                                 neighbor_level,
                                 /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          failures += DD::check_for_broken_coupling_intersections(coupling_glue);
        }
      }
    }
    return failures;
  } // ... count_wrong_intersections(...)

  std::unique_ptr<GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<DD::Glued<MacroGridType, LocalGridType, local_layer>> dd_grid_;
}; // struct GluedDdGridTest


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_GRID_GLUE
#endif // DUNE_XT_GRID_TEST_DD_GLUED_HH
