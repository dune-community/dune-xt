// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2020)
//   René Fritze     (2015 - 2020)
//   Tobias Leibner  (2016 - 2017, 2019 - 2020)

#include "config.h"

#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <vector>
#include <map>
#include <string>
#include <sstream>

#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/test/common.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>


void busywait(const size_t ms)
{
  // "round" up to next full 10 ms to align with native timer res
  const size_t milliseconds = (ms / 10) * 10 + 10;
  timeval start, end;
  gettimeofday(&start, NULL);
  do {
    gettimeofday(&end, NULL);
  } while (((end.tv_sec - start.tv_sec) * 1e6) + ((end.tv_usec - start.tv_usec)) < milliseconds * 1000);
} // ... busywait(...)


Dune::XT::Common::Configuration DXTC_TEST_CONFIG_SUB(const std::string& sub_key)
{
  const auto key = Dune::XT::Common::Test::get_unique_test_name() + "." + sub_key;
  EXPECT_TRUE(DXTC_CONFIG.has_sub(key));
  return DXTC_CONFIG.sub(key);
}


namespace Dune {
namespace XT {
namespace Common {
namespace Test {


std::string get_unique_test_name()
{
  const ::testing::TestInfo* const test_info = ::testing::UnitTest::GetInstance()->current_test_info();
  DUNE_THROW_IF(!test_info, Exceptions::you_are_using_this_wrong, "You may only use this in the context of gtest!");
  // For grids, put the names used in <dune/xt/grid/grids.hh> here.
  static std::map<std::string, std::string> replacements{
      {"Dune::ALUGrid<2, 2, (Dune::ALUGridElementType)0, (Dune::ALUGridRefinementType)0, Dune::ALUGridMPIComm>",
       "ALU_2D_SIMPLEX_CONFORMING"},
      {"Dune::ALUGrid<2, 2, (Dune::ALUGridElementType)0, (Dune::ALUGridRefinementType)1, Dune::ALUGridMPIComm>",
       "ALU_2D_SIMPLEX_NONCONFORMING"},
      {"Dune::ALUGrid<2, 2, (Dune::ALUGridElementType)1, (Dune::ALUGridRefinementType)1, Dune::ALUGridMPIComm>",
       "ALU_2D_CUBE"},
      {"Dune::ALUGrid<3, 3, (Dune::ALUGridElementType)0, (Dune::ALUGridRefinementType)0, Dune::ALUGridMPIComm>",
       "ALU_3D_SIMPLEX_CONFORMING"},
      {"Dune::ALUGrid<3, 3, (Dune::ALUGridElementType)0, (Dune::ALUGridRefinementType)1, Dune::ALUGridMPIComm>",
       "ALU_3D_SIMPLEX_NONCONFORMING"},
      {"Dune::ALUGrid<3, 3, (Dune::ALUGridElementType)1, (Dune::ALUGridRefinementType)1, Dune::ALUGridMPIComm>",
       "ALU_3D_CUBE"},
      {"Dune::OneDGrid", "ONED_1D"},
      {"Dune::UGGrid<2>", "UG_2D"},
      {"Dune::UGGrid<3>", "UG_3D"},
      {"Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1> >", "YASP_1D_EQUIDISTANT_OFFSET"},
      {"Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2> >", "YASP_2D_EQUIDISTANT_OFFSET"},
      {"Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3> >", "YASP_3D_EQUIDISTANT_OFFSET"},
      {"Dune::YaspGrid<4, Dune::EquidistantOffsetCoordinates<double, 4> >", "YASP_4D_EQUIDISTANT_OFFSET"}};
  auto replace_if = [&](const auto& id) {
    if (replacements.count(std::string(id)) > 0)
      return replacements.at(std::string(id));
    else
      return std::string(id);
  };
  std::string result;
  const auto* test_case_name = test_info->test_case_name();
  if (test_case_name) {
    result += replace_if(test_case_name);
  }
  const auto* name = test_info->name();
  if (name) {
    result += result.empty() ? "" : "." + replace_if(name);
  }
  const auto* type_param = test_info->type_param();
  if (type_param) {
    result += result.empty() ? "" : "." + replace_if(type_param);
  }
  const auto* value_param = test_info->value_param();
  if (value_param) {
    result += result.empty() ? "" : "." + replace_if(value_param);
  }
  std::replace(result.begin(), result.end(), '/', '.');
  return result;
}


} // namespace Test
} // namespace Common
namespace Test {
namespace internal {


std::pair<size_t, ssize_t> convert_to_scientific(const double number, const size_t precision)
{
  // see http://www.mathworks.com/matlabcentral/newsreader/view_thread/151859
  const double sign = (number > 0.0) ? 1.0 : -1.0;
  const double tmp = std::log10(std::abs(number));
  ssize_t exponent = ((tmp > 0) ? 1 : -1) * std::floor(std::abs(tmp));
  double coefficient = sign * std::pow(10.0, tmp - exponent);
  if (std::abs(coefficient) < 1.0) {
    coefficient *= 10.0;
    exponent -= 1;
  }
  const double factor = std::pow(10, precision);
  return std::make_pair(size_t(std::round(factor * coefficient)), exponent);
} // ... convert_to_scientific(...)


} // namespace internal


void check_eoc_study_for_success(
    const Common::Configuration& expected_results,
    const std::map<std::string, std::map<std::string, std::map<size_t, double>>>& actual_results,
    const double& zero_tolerance)
{
  auto get_keys = [](const auto& map) {
    std::vector<std::string> result;
    for (const auto& ele : map)
      result.push_back(Common::to_string(ele.first));
    return result;
  };
  auto clean_up = [](const auto& str) {
    std::string cleaned(str);
    std::replace(cleaned.begin(), cleaned.end(), '/', '_');
    std::replace(cleaned.begin(), cleaned.end(), ' ', '_');
    return cleaned;
  };
  auto get_matching_key = [&clean_up](const auto& map, const auto& str) {
    size_t matches = 0;
    std::string matching_key;
    for (const auto& entry : map) {
      if (clean_up(entry.first) == str) {
        matching_key = entry.first;
        ++matches;
      }
    }
    EXPECT_EQ(1, matches) << "There were multiple keys in the actual results which (after clean up) matched a key in "
                             "the expected results, this must not be!";
    return matching_key;
  };
  EXPECT_GT(expected_results.getSubKeys().size(), 0)
      << "you have to provide at least one category!"
      << "\nThe actual results for instance contain the categories '" << get_keys(actual_results) << "'."
      << "\n\n"
      << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
  for (const auto& category : expected_results.getSubKeys()) {
    EXPECT_GT(actual_results.count(category), 0) << "missing data for category " << category << "!"
                                                 << "\n\n"
                                                 << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
    const auto expected_category_data = expected_results.sub(category);
    const auto& actual_category_data = actual_results.at(get_matching_key(actual_results, category));
    EXPECT_GT(expected_category_data.getValueKeys().size(), 0)
        << "you have to provide at least one type for category " << category
        << "!\nThe actual_results for instance contain the types '" << get_keys(actual_category_data)
        << "' for this category."
        << "\n\n"
        << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
    for (const auto& type : expected_category_data.getValueKeys()) {
      EXPECT_GT(actual_category_data.count(get_matching_key(actual_category_data, type)), 0)
          << "missing data for category " << category << " and type " << type << "!"
          << "\n\n"
          << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
      const auto expected_values = expected_category_data.template get<std::vector<double>>(type);
      const auto& actual_type_data = actual_category_data.at(get_matching_key(actual_category_data, type));
      for (size_t level = 0; level < expected_values.size(); ++level) {
        EXPECT_GT(actual_type_data.count(level), 0)
            << "missing data for category " << category << ", type " << type << " and level " << level << "!"
            << "\n\n"
            << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
        const auto& actual_value = actual_type_data.at(level);
        const auto& expected_value = expected_values[level];
        if (std::abs(expected_value) < zero_tolerance) {
          const auto abs_of_actual_value = std::abs(actual_value);
          EXPECT_LT(abs_of_actual_value, zero_tolerance)
              << "          Expected result is interpreted as zero and result is not close enough to zero!\n"
              << "          level:    " << level << "\n"
              << "          category: " << category << "\n"
              << "          type:     " << type << "\n"
              << "          zero_tolerance: " << zero_tolerance << "\n"
              << "          actual_value   = " << actual_value << "\n"
              << "          expected_value = " << expected_value << "\n\n"
              << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
        } else {
          const auto actual_result = internal::convert_to_scientific(actual_value, 2);
          const auto expected_result = internal::convert_to_scientific(expected_value, 2);
          const auto actual_exponent = actual_result.second;
          const auto expected_exponent = expected_result.second;
          EXPECT_EQ(expected_exponent, actual_exponent)
              << "          Exponent comparison (in scientific notation, precision 2) failed for\n"
              << "          level:    " << level << "\n"
              << "          category: " << category << "\n"
              << "          type:     " << type << "\n"
              << "          actual_value   = " << actual_value << "\n"
              << "          expected_value = " << expected_value << "\n\n"
              << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
          if (actual_exponent != expected_exponent)
            return;
          const auto actual_coefficient = actual_result.first;
          const auto expected_coefficient = expected_result.first;
          EXPECT_EQ(expected_coefficient, actual_coefficient)
              << "          Coefficient comparison (in scientific notation, precision 2) failed for\n"
              << "          level:    " << level << "\n"
              << "          category: " << category << "\n"
              << "          type:     " << type << "\n"
              << "          actual_value   = " << actual_value << "\n"
              << "          expected_value = " << expected_value << "\n\n"
              << "get_unique_test_name() = " << Common::Test::get_unique_test_name();
        }
      }
    }
  }
} // ... check_eoc_study_for_success(...)


void print_collected_eoc_study_results(const std::map<std::string, std::vector<double>>& results, std::ostream& out)
{
  if (results.empty())
    DUNE_THROW(Common::Exceptions::wrong_input_given, "Given results must not be empty!");
  std::vector<std::string> actually_used_norms;
  for (const auto& element : results)
    actually_used_norms.push_back(element.first);
  out << "if (type == \"" << actually_used_norms[0] << "\")\n";
  out << "  return " << results.at(actually_used_norms[0]) << ";\n";
  for (size_t ii = 1; ii < actually_used_norms.size(); ++ii) {
    out << "else if (type == \"" << actually_used_norms[ii] << "\")\n";
    out << "  return " << results.at(actually_used_norms[ii]) << ";\n";
  }
  out << "else\n";
  out << "  EXPECT_TRUE(false) << \"test results missing for type: \" << type;\n";
  out << "return {};" << std::endl;
} // ... print_collected_eoc_study_results(...)


unsigned int grid_elements()
{
  return DXTC_CONFIG.has_key("test.gridelements") // <- doing this so complicated to
             ? DXTC_CONFIG.get<unsigned int>(
                 "test.gridelements", 3u, Common::ValidateLess<unsigned int>(2u)) //    silence the WARNING: ...
             : 3u;
} // ... grid_elements(...)


} // namespace Test
} // namespace XT
} // namespace Dune
