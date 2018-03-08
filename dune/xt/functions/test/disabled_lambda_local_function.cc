// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#include <dune/xt/common/test/main.hxx>

GTEST_TEST(DISABLED_LocalLambdaFluxFunction, function_is_not_yet_up_to_date_with_new_interface)
{
}

#if 0
struct LocalLambdaFluxFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  void check() const
  {
    typedef TESTFUNCTIONTYPE F;

    F f(
        [](const typename F::EntityType& /*entity*/,
           const typename F::DomainType& xx,
           const XT::Common::Parameter& mu) {
          typename F::RangeType ret(std::pow(xx[0], mu.get("power").at(0)));
          return ret;
        },
        /*order=*/8, // <- dos not make much sense in this setting
        XT::Common::ParameterType("power", 1),
        "x_power_p");

    auto grid = XT::Grid::make_cube_grid<GRIDTYPE>();

    for (auto&& entity : elements(grid.leaf_view())) {
      auto xx_global = entity.geometry().center();
      auto xx_local = entity.geometry().local(xx_global);
      auto f_val = f.local_function(entity)->evaluate(xx_local, {"power", 2.});
      ASSERT_EQ(std::pow(xx_local[0], 2.), get_first_entry(f_val));
    }
  }
};
#endif // 0
