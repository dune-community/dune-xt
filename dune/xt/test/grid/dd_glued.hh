// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_GRID_TEST_DD_GLUED_HH
#define DUNE_XT_GRID_TEST_DD_GLUED_HH
#if HAVE_DUNE_GRID_GLUE

#  include <map>
#  include <set>
#  include <sstream>
#  include <vector>

#  include <dune/grid/common/rangegenerators.hh>

#  include <dune/xt/common/memory.hh>
#  include <dune/xt/test/gtest/gtest.h>
#  include <dune/xt/grid/grids.hh>
#  include <dune/xt/grid/gridprovider.hh>
#  include <dune/xt/grid/dd/glued.hh>


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


template <class M, class L>
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
    static constexpr Layers type = Layers::level;
  };

#  if HAVE_ALBERTA

  template <int d, int dw, bool anything>
  struct get_local_layer<AlbertaGrid<d, dw>, anything>
  {
    static constexpr Layers type = Layers::leaf;
  };

#  endif

  static constexpr Layers local_layer = get_local_layer<LocalGridType>::type;

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

  std::unique_ptr<GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<DD::Glued<MacroGridType, LocalGridType, local_layer>> dd_grid_;
}; // struct GluedDdGridTest


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_GRID_GLUE
#endif // DUNE_XT_GRID_TEST_DD_GLUED_HH
