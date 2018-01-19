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

#ifndef DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif

#include <dune/xt/common/test/main.hxx>

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/flattop.hh>

#include "functions.hh"

using namespace Dune;
using namespace Dune::XT;

template <class G>
class FlatTopFunctionType
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const size_t d = G::dimension;
  typedef double R;
  static const size_t r = 1;
  static const size_t rC = 1;

public:
  typedef Functions::FlatTopFunction<E, D, d, R, r, rC> value;
}; // struct FlatTopFunctionType

template <class DimDomain>
class FlatTopFunctionTest
    : public FunctionTest<
          typename FlatTopFunctionType<YaspGrid<DimDomain::value,
                                                EquidistantOffsetCoordinates<double, DimDomain::value>>>::value>
{
protected:
  typedef YaspGrid<DimDomain::value, EquidistantOffsetCoordinates<double, DimDomain::value>> GridType;
  typedef typename FlatTopFunctionType<GridType>::value FunctionType;

  static std::shared_ptr<GridType> create_grid()
  {
    return XT::Grid::make_cube_grid<GridType>(0, 3, 12).grid_ptr();
  }

  template <class P, class V, class L, class R, class D, class E>
  static void check(const P& point, const V& value, const L& left, const R& right, const D& delta, const E& top_value)
  {
    if (Common::FloatCmp::lt(point, left - delta) || Common::FloatCmp::gt(point, right + delta)) {
      // outside
      EXPECT_EQ(0.0, value) << point;
    } else if (Common::FloatCmp::ge(point, left + delta) && Common::FloatCmp::le(point, right - delta)) {
      // inside
      EXPECT_EQ(top_value, value) << point;
    } else {
      // boundary layer
      if (top_value > 0.0) {
        EXPECT_GE(top_value, value) << point[0];
        EXPECT_LE(0.0, value) << point[0];
      } else {
        EXPECT_LE(top_value, value) << point[0];
        EXPECT_GE(0.0, value) << point[0];
      }
    }
  } // ... check(...)
}; // class FlatTopFunctionTest

typedef testing::Types<Int<1>, Int<2>, Int<3>> DimDomains;

TYPED_TEST_CASE(FlatTopFunctionTest, DimDomains);
TYPED_TEST(FlatTopFunctionTest, static_interface_check)
{
  this->static_interface_check();
}
TYPED_TEST(FlatTopFunctionTest, static_create_check)
{
  this->static_create_check();
}
TYPED_TEST(FlatTopFunctionTest, dynamic_interface_check)
{
  this->dynamic_interface_check(*(TestFixture::FunctionType::create()), *(this->create_grid()));
}
TYPED_TEST(FlatTopFunctionTest, evaluate_check)
{
  auto grid_ptr = this->create_grid();
  typedef FieldVector<double, TypeParam::value> DomainType;
  const DomainType left(1);
  const DomainType right(2);
  const DomainType delta(1e-6);
  const double value = 20;
  typename TestFixture::FunctionType func(left, right, delta, value, "bar");
  func.visualize(grid_ptr->leafGridView(), "dim_" + Common::to_string(int(TypeParam::value)));
  for (auto&& entity : elements(grid_ptr->leafGridView())) {
    const auto local_func = func.local_function(entity);
    const auto& quadrature = QuadratureRules<double, TypeParam::value>::rule(
        entity.type(), boost::numeric_cast<int>(local_func->order() + 2));
    for (const auto& element : quadrature) {
      const auto& local_point = element.position();
      const auto point = entity.geometry().global(local_point);
      const auto val = local_func->evaluate(local_point);
      this->check(point, val, left, right, delta, value);
    }
  }
}
