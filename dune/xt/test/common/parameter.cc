// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#include <dune/xt/common/test/main.hxx>

#include <iostream>

#include <dune/xt/common/parameter.hh>

using namespace Dune::XT::Common;


GTEST_TEST(ParameterType, creation_and_report_and_ostreamout)
{
  std::map<ParameterType, std::string> expected_values = {
      {ParameterType(), "ParameterType({})"},
      {ParameterType("single_key"), "ParameterType({single_key: 1})"},
      {ParameterType("single_key_with_size_17", 17), "ParameterType({single_key_with_size_17: 17})"},
      {ParameterType({{"first_key_with_size_1", 1}, {"second_key_with_size_2", 2}}),
       "ParameterType({first_key_with_size_1: 1,\n               second_key_with_size_2: 2})"}};

  for (const auto& element : expected_values) {
    EXPECT_EQ(element.second, element.first.report());
    std::stringstream ss;
    ss << element.first;
    EXPECT_EQ(element.second, ss.str());
  }
} // ParameterType, creation_and_report_and_ostreamout


GTEST_TEST(ParameterType, equality_comparison)
{
  using P = std::pair<ParameterType, std::vector<ParameterType>>;
  std::list<P> those_are_equal;
  those_are_equal.emplace_back(P(ParameterType(), {ParameterType()}));
  those_are_equal.emplace_back(P(ParameterType("__unspecified__"),
                                 {ParameterType("__unspecified__"),
                                  ParameterType("__unspecified__", 1),
                                  ParameterType("And", 1),
                                  ParameterType("Any", 1),
                                  ParameterType("other", 1),
                                  ParameterType("singLe", 1),
                                  ParameterType("word", 1)}));
  those_are_equal.emplace_back(P(ParameterType("__unspecified__", 3),
                                 {ParameterType("__unspecified__", 3),
                                  ParameterType("And", 3),
                                  ParameterType("Any", 3),
                                  ParameterType("other", 3),
                                  ParameterType("singLe", 3),
                                  ParameterType("word", 3)}));
  those_are_equal.emplace_back(P(ParameterType("scalar"),
                                 {ParameterType("scalar"),
                                  ParameterType("scalar", 1),
                                  ParameterType("__unspecified__"),
                                  ParameterType("__unspecified__", 1)}));
  those_are_equal.emplace_back(
      P(ParameterType("vector", 17), {ParameterType("vector", 17), ParameterType("__unspecified__", 17)}));
  for (const auto& element : those_are_equal) {
    const auto& target = element.first;
    for (const auto& source : element.second)
      EXPECT_EQ(target, source);
  }
} // ParameterType, equality_comparison


GTEST_TEST(ParameterType, inequality_comparison)
{
  using P = std::pair<ParameterType, std::vector<ParameterType>>;
  std::list<P> those_are_not_equal;
  those_are_not_equal.emplace_back(P(ParameterType(),
                                     {ParameterType("__unspecified__"),
                                      ParameterType("scalar"),
                                      ParameterType("scalar", 1),
                                      ParameterType("vector", 17)}));
  those_are_not_equal.emplace_back(
      P(ParameterType("__unspecified__"),
        {ParameterType(), ParameterType("__unspecified__", 3), ParameterType("vector", 17)}));
  those_are_not_equal.emplace_back(P(ParameterType("__unspecified__", 3),
                                     {ParameterType(),
                                      ParameterType("__unspecified__"),
                                      ParameterType("scalar"),
                                      ParameterType("scalar", 1),
                                      ParameterType("vector", 17)}));
  those_are_not_equal.emplace_back(P(ParameterType("scalar"),
                                     {ParameterType(),
                                      ParameterType("another_scalar"),
                                      ParameterType("scalar", 17),
                                      ParameterType("__unspecified__", 3),
                                      ParameterType("vector", 17)}));
  those_are_not_equal.emplace_back(P(ParameterType("vector", 17),
                                     {ParameterType(),
                                      ParameterType("__unspecified__"),
                                      ParameterType("__unspecified__", 3),
                                      ParameterType("scalar"),
                                      ParameterType("another_vector", 17),
                                      ParameterType("vector", 3)}));
  for (const auto& element : those_are_not_equal) {
    const auto& target = element.first;
    for (const auto& source : element.second)
      EXPECT_NE(target, source);
  }
} // ParameterType, inequality_comparison


GTEST_TEST(Parameter, creation_and_report_and_ostreamout)
{
  std::map<Parameter, std::string> expected_values = {
      {Parameter(), "Parameter({})"},
      {Parameter(1.), "Parameter({__unspecified__: 1})"},
      {Parameter({1., 2.}), "Parameter({__unspecified__: [1 2]})"},
      {Parameter("scalar_value", 17.), "Parameter({scalar_value: 17})"},
      {Parameter("vector_value", {17., 42.}), "Parameter({vector_value: [17 42]})"},
      {Parameter({{"first_scalar_value", {1.}}, {"second_vector_value", {1., 2.}}}),
       "Parameter({first_scalar_value: 1,\n           second_vector_value: [1 2]})"}};
  for (const auto& element : expected_values) {
    EXPECT_EQ(element.second, element.first.report());
    std::stringstream ss;
    ss << element.first;
    EXPECT_EQ(element.second, ss.str());
  }
} // Parameter, creation_and_report_and_ostreamout
