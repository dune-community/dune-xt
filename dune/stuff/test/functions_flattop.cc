// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif

#include "main.hxx"

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/stuff/common/disable_warnings.hh>
#if HAVE_DUNE_GRID
#include <dune/grid/sgrid.hh>
#endif
#include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/playground/functions/flattop.hh>

#include "function.hh"

#if HAVE_DUNE_GRID

using namespace Dune;
using namespace Stuff;


template <class G>
class FlatTopFunctionType
{
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const unsigned int d = G::dimension;
  typedef double R;
  static const unsigned int r  = 1;
  static const unsigned int rC = 1;

public:
  typedef Functions::FlatTop<E, D, d, R, r, rC> value;
}; // struct FlatTopFunctionType


template <class DimDomain>
class FlatTopFunctionTest
    : public FunctionTest<typename FlatTopFunctionType<SGrid<DimDomain::value, DimDomain::value>>::value>
{
protected:
  typedef SGrid<DimDomain::value, DimDomain::value> GridType;
  typedef typename FlatTopFunctionType<GridType>::value FunctionType;

  static std::shared_ptr<GridType> create_grid()
  {
    return Stuff::Grid::Providers::Cube<GridType>(0.0, 3.0, 12).grid_ptr();
  }

  template <class... Args>
  static std::unique_ptr<FunctionType> create(Args&&... args)
  {
    return std::unique_ptr<FunctionType>(new FunctionType(std::forward<Args>(args)...));
  }

  template <class P, class V, class L, class R, class D, class E>
  static void check(const P& point, const V& value, const L& left, const R& right, const D& delta, const E& top_value)
  {
    if (Common::FloatCmp::lt(point, left - delta) || Common::FloatCmp::gt(point, right + delta)) {
      // outside
      EXPECT_EQ(0.0, value) << point;
    } else if (Common::FloatCmp::ge(point, left) && Common::FloatCmp::le(point, right)) {
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
TYPED_TEST(FlatTopFunctionTest, dynamic_interface_check)
{
  this->dynamic_interface_check(*(this->create()), *(this->create_grid()));
}
TYPED_TEST(FlatTopFunctionTest, copy_check)
{
  this->copy_check(*(this->create()));
}
TYPED_TEST(FlatTopFunctionTest, evaluate_check)
{
  auto grid_ptr = this->create_grid();
  typedef FieldVector<double, TypeParam::value> DomainType;
  const DomainType left(1.0);
  const DomainType right(2.0);
  const DomainType delta(0.25);
  const double value = 1.75;
  auto func          = this->create(left, right, delta, value, "bar");
  //  func->visualize(grid_ptr->leafGridView(), "foo");
  for (const auto& entity : Stuff::Common::viewRange(grid_ptr->leafGridView())) {
    const auto local_func  = func->local_function(entity);
    const auto& quadrature = QuadratureRules<double, TypeParam::value>::rule(
        entity.type(), boost::numeric_cast<int>(local_func->order() + 2));
    for (const auto& element : quadrature) {
      const auto& local_point = element.position();
      const auto point        = entity.geometry().global(local_point);
      const auto val = local_func->evaluate(local_point);
      this->check(point, val, left, right, delta, value);
    }
  }
}


#else // HAVE_DUNE_GRID


TYPED_TEST(DISABLED_FlatTopFunctionTest, static_interface_check){};


#endif // HAVE_DUNE_GRID
