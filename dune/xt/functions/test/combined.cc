// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2015 - 2016, 2018)
//   Tobias Leibner  (2016)

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/base/combined.hh>
#include <dune/xt/functions/constant.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

template <class G>
class DifferenceFunctionType
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef double R;
  static const size_t r = 1;
  static const size_t rC = 1;

public:
  typedef Functions::ConstantFunction<E, D, d, R, r, rC> ConstantFunctionType;
  typedef Functions::DifferenceFunction<ConstantFunctionType, ConstantFunctionType> value;
}; // struct DifferenceFunctionType

template <class DimDomain>
class DifferenceFunctionTest
    : public FunctionTest<
          typename DifferenceFunctionType<YaspGrid<DimDomain::value,
                                                   EquidistantOffsetCoordinates<double, DimDomain::value>>>::value>
{
protected:
  typedef YaspGrid<DimDomain::value, EquidistantOffsetCoordinates<double, DimDomain::value>> GridType;
  typedef typename DifferenceFunctionType<GridType>::value FunctionType;

  static std::shared_ptr<GridType> create_grid()
  {
    return XT::Grid::make_cube_grid<GridType>(0.0, 1.0, 4).grid_ptr();
  }

  static std::unique_ptr<FunctionType> create(const double ll, const double rr)
  {
    typedef typename DifferenceFunctionType<YaspGrid<DimDomain::value,
                                                     EquidistantOffsetCoordinates<double, DimDomain::value>>>::
        ConstantFunctionType ConstantFunctionType;
    auto left = std::make_shared<ConstantFunctionType>(ll);
    auto right = std::make_shared<ConstantFunctionType>(rr);
    return std::unique_ptr<FunctionType>(new FunctionType(left, right));
  } // ... create(...)
}; // class DifferenceFunctionTest

typedef testing::Types<Int<1>, Int<2>, Int<3>> DimDomains;

TYPED_TEST_CASE(DifferenceFunctionTest, DimDomains);
TYPED_TEST(DifferenceFunctionTest, static_interface_check)
{
  this->static_interface_check();
}
TYPED_TEST(DifferenceFunctionTest, dynamic_interface_check)
{
  this->dynamic_interface_check(*(this->create(1.0, 1.0)), *(this->create_grid()));
}
TYPED_TEST(DifferenceFunctionTest, evaluate_check)
{
  auto grid_ptr = this->create_grid();
  auto func = this->create(1.0, 2.0);
  //  func->visualize(grid_ptr->leafGridView(), "foo");
  for (auto&& entity : elements(grid_ptr->leafGridView())) {
    const auto local_func = func->local_function(entity);
    const auto& quadrature = QuadratureRules<double, TypeParam::value>::rule(
        entity.type(), boost::numeric_cast<int>(local_func->order() + 2));
    for (const auto& element : quadrature) {
      const auto& local_point = element.position();
      const auto val = local_func->evaluate(local_point);
      EXPECT_EQ(val[0], -1.0);
    }
  }
} // DifferenceFunctionTest, evaluate_check
