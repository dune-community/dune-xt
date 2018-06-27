// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014 - 2017)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <dune/xt/common/exceptions.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#if HAVE_DUNE_XT_LA
#include <dune/xt/la/container/eye-matrix.hh>
#endif // HAVE_DUNE_XT_LA

#include <dune/xt/functions/affine.hh>
#include <dune/xt/functions/test/functions.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::Functions;

#if HAVE_DUNE_XT_LA
struct AffineFunctionTest : public FunctionTest<TESTFUNCTIONTYPE>
{
  void check()
  {
    using MatrixType = typename TESTFUNCTIONTYPE::FieldMatrixType;
    static constexpr size_t dimDomain = TESTFUNCTIONTYPE::dimDomain;
    static constexpr size_t dimRange = TESTFUNCTIONTYPE::dimRange;
    static constexpr size_t dimRangeCols = TESTFUNCTIONTYPE::dimRangeCols;
    auto grid = XT::Grid::make_cube_grid<GRIDTYPE>();
    auto unit_matrix = XT::LA::eye_matrix<MatrixType>(dimRange, dimDomain);
    auto matrices = FieldVector<MatrixType, dimRangeCols>(unit_matrix);
    for (size_t rC = 0; rC < dimRangeCols; ++rC)
      matrices[rC] *= rC + 1;
    const auto testfunction = TESTFUNCTIONTYPE(matrices);
    for (auto&& entity : elements(grid.leaf_view())) {
      auto xx_global = entity.geometry().center();
      auto xx_local = entity.geometry().local(xx_global);
      FieldVector<MatrixType, dimRangeCols> jacobian;
      auto ret = testfunction.local_function(entity)->evaluate(xx_local);
      for (size_t ii = 0; ii < std::min(dimRange, dimDomain); ++ii) {
        FieldVector<double, dimRangeCols> ret_row(ret[ii]);
        for (size_t rC = 0; rC < dimRangeCols; ++rC)
          EXPECT_DOUBLE_EQ(xx_global[ii] * (rC + 1), ret_row[rC]);
      }
      for (size_t ii = std::min(dimRange, dimDomain) + 1; ii < dimRange; ++ii) {
        FieldVector<double, dimRangeCols> ret_row(ret[ii]);
        for (size_t rC = 0; rC < dimRangeCols; ++rC)
          EXPECT_DOUBLE_EQ(0., ret_row[rC]);
      }
      testfunction.local_function(entity)->jacobian(xx_local, jacobian);
      for (size_t rC = 0; rC < dimRangeCols; ++rC)
        ASSERT_TRUE(XT::Common::FloatCmp::eq(jacobian[rC], matrices[rC]));
    }
  }
};

TEST_F(AffineFunctionTest, creation_and_evalution)
{
  this->check();
}
#endif // HAVE_DUNE_XT_LA
