// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)


#include <dune/xt/test/main.hxx> // <- has to come first, include config.h!
#include <dune/xt/functions/grid-function.hh>

using namespace Dune::XT;

using E1 = Grid::extract_entity_t<YASP_1D_EQUIDISTANT_OFFSET>;
using E2 = Grid::extract_entity_t<YASP_2D_EQUIDISTANT_OFFSET>;

#define VAL(r, rC) Functions::RangeTypeSelector<double, r, rC>::return_type(0)


template <class L, class R>
void test_difference(const L& left, const R& right)
{
  auto ll = left.local_function();
  auto lr = right.local_function();
  (*ll) - (*lr);
}


template <class L, class R>
void test_fraction(const L& left, const R& right)
{
  auto ll = left.local_function();
  auto lr = right.local_function();
  (*ll) / (*lr);
}


template <class L, class R>
void test_product(const L& left, const R& right)
{
  auto ll = left.local_function();
  auto lr = right.local_function();
  (*ll) * (*lr);
}


template <class L, class R>
void test_sum(const L& left, const R& right)
{
  auto ll = left.local_function();
  auto lr = right.local_function();
  (*ll) + (*lr);
}


GTEST_TEST(grid_function_interface_operators, main)
{
  test_difference(Functions::GridFunction<E1, 1, 1>(VAL(1, 1)), Functions::GridFunction<E1, 1, 1>(VAL(1, 1)));
  test_difference(Functions::GridFunction<E2, 1, 1>(VAL(1, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_difference(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 1, 2>(VAL(1, 2)));
  test_difference(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 1, 3>(VAL(1, 3)));
  test_difference(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));
  test_difference(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 2, 2>(VAL(2, 2)));
  test_difference(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 2, 3>(VAL(2, 3)));
  test_difference(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));
  test_difference(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 3, 2>(VAL(3, 2)));
  test_difference(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 3, 3>(VAL(3, 3)));

  test_sum(Functions::GridFunction<E1, 1, 1>(VAL(1, 1)), Functions::GridFunction<E1, 1, 1>(VAL(1, 1)));
  test_sum(Functions::GridFunction<E2, 1, 1>(VAL(1, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_sum(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 1, 2>(VAL(1, 2)));
  test_sum(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 1, 3>(VAL(1, 3)));
  test_sum(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));
  test_sum(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 2, 2>(VAL(2, 2)));
  test_sum(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 2, 3>(VAL(2, 3)));
  test_sum(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));
  test_sum(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 3, 2>(VAL(3, 2)));
  test_sum(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 3, 3>(VAL(3, 3)));

  test_fraction(Functions::GridFunction<E1, 1, 1>(VAL(1, 1)), Functions::GridFunction<E1, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 1, 1>(VAL(1, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_fraction(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));

  test_product(Functions::GridFunction<E1, 1, 1>(VAL(1, 1)), Functions::GridFunction<E1, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 1, 1>(VAL(1, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));
  test_product(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 1, 1>(VAL(1, 1)));

  test_product(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));
  test_product(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 2, 2>(VAL(2, 2)));
  test_product(Functions::GridFunction<E2, 1, 2>(VAL(1, 2)), Functions::GridFunction<E2, 2, 3>(VAL(2, 3)));

  test_product(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));
  test_product(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 3, 2>(VAL(3, 2)));
  test_product(Functions::GridFunction<E2, 1, 3>(VAL(1, 3)), Functions::GridFunction<E2, 3, 3>(VAL(3, 3)));

  test_product(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 1, 2>(VAL(1, 2)));
  test_product(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 1, 3>(VAL(1, 3)));
  test_product(Functions::GridFunction<E2, 2, 1>(VAL(2, 1)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));

  test_product(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));
  test_product(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 2, 2>(VAL(2, 2)));
  test_product(Functions::GridFunction<E2, 2, 2>(VAL(2, 2)), Functions::GridFunction<E2, 2, 3>(VAL(2, 3)));

  test_product(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));
  test_product(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 3, 2>(VAL(3, 2)));
  test_product(Functions::GridFunction<E2, 2, 3>(VAL(2, 3)), Functions::GridFunction<E2, 3, 3>(VAL(3, 3)));

  test_product(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 1, 2>(VAL(1, 2)));
  test_product(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 1, 3>(VAL(1, 3)));
  test_product(Functions::GridFunction<E2, 3, 1>(VAL(3, 1)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));

  test_product(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 2, 1>(VAL(2, 1)));
  test_product(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 2, 2>(VAL(2, 2)));
  test_product(Functions::GridFunction<E2, 3, 2>(VAL(3, 2)), Functions::GridFunction<E2, 2, 3>(VAL(2, 3)));

  test_product(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 3, 1>(VAL(3, 1)));
  test_product(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 3, 2>(VAL(3, 2)));
  test_product(Functions::GridFunction<E2, 3, 3>(VAL(3, 3)), Functions::GridFunction<E2, 3, 3>(VAL(3, 3)));
}
