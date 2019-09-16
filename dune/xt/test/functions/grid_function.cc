// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)

#include <dune/xt/common/test/main.hxx> // <- Has to come first, includes the config.h!

#include <dune/xt/common/unused.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>
#include <dune/xt/functions/grid-function.hh>

using namespace Dune;
using namespace Dune::XT::Functions;

using G = CUBEGRID_2D;
using E = XT::Grid::extract_entity_t<G>;
static const constexpr size_t d = G::dimension;
static const constexpr size_t r = 2;
static const constexpr size_t rC = 3;

// from here on, the code should work for any E and d to allow for grid parametrization


std::string accepts_scalar_grid_function(GridFunction<E> func)
{
  return func.name();
}

std::string accepts_vector_grid_function(GridFunction<E, r> func)
{
  return func.name();
}

std::string accepts_square_matrix_grid_function(GridFunction<E, r, r> func)
{
  return func.name();
}

std::string accepts_matrix_grid_function(GridFunction<E, r, rC> func)
{
  return func.name();
}


// ScalarGridFunction

GTEST_TEST(ScalarGridFunction, constructible_from_double)
{
  GridFunction<E> DXTC_UNUSED(func){1.};
}
GTEST_TEST(ScalarGridFunction, convertible_from_double)
{
  accepts_scalar_grid_function(1.);
}
GTEST_TEST(ScalarGridFunction, constructible_from_FieldVector)
{
  GridFunction<E> DXTC_UNUSED(func){FieldVector<double, 1>()};
}
GTEST_TEST(ScalarGridFunction, convertible_from_FieldVector)
{
  accepts_scalar_grid_function(FieldVector<double, 1>());
}
GTEST_TEST(ScalarGridFunction, constructible_from_XT_Common_FieldVector)
{
  GridFunction<E> DXTC_UNUSED(func){XT::Common::FieldVector<double, 1>()};
}
GTEST_TEST(ScalarGridFunction, convertible_from_XT_Common_FieldVector)
{
  accepts_scalar_grid_function(XT::Common::FieldVector<double, 1>());
}
GTEST_TEST(ScalarGridFunction, constructible_from_function_rvalueref)
{
  const GenericFunction<d> scalar_function{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  GridFunction<E> func{scalar_function};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(ScalarGridFunction, convertible_from_function_rvalueref)
{
  const GenericFunction<d> scalar_function{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  auto nm = accepts_scalar_grid_function(scalar_function);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(ScalarGridFunction, constructible_from_function_lvalueptr)
{
  GridFunction<E> func{new GenericFunction<d>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(ScalarGridFunction, convertible_from_function_lvalueptr)
{
  accepts_scalar_grid_function(new GenericFunction<d>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"});
}
GTEST_TEST(ScalarGridFunction, constructible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E> scalar_grid_function{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  GridFunction<E> func{scalar_grid_function};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(ScalarGridFunction, convertible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E> scalar_grid_function{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  auto nm = accepts_scalar_grid_function(scalar_grid_function);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(ScalarGridFunction, constructible_from_grid_function_lvalueptr)
{
  GridFunction<E> func{new GenericGridFunction<E>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(ScalarGridFunction, convertible_from_grid_function_lvalueptr)
{
  accepts_scalar_grid_function(
      new GenericGridFunction<E>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"});
}


// VectorGridFunction

GTEST_TEST(VectorGridFunction, constructible_from_FieldVector)
{
  GridFunction<E, r> DXTC_UNUSED(func){FieldVector<double, r>()};
}
GTEST_TEST(VectorGridFunction, convertible_from_FieldVector)
{
  accepts_vector_grid_function(FieldVector<double, r>());
}
GTEST_TEST(VectorGridFunction, constructible_from_XT_Common_FieldVector)
{
  GridFunction<E, r> DXTC_UNUSED(func){XT::Common::FieldVector<double, r>()};
}
GTEST_TEST(VectorGridFunction, convertible_from_XT_Common_FieldVector)
{
  accepts_vector_grid_function(XT::Common::FieldVector<double, r>());
}
GTEST_TEST(VectorGridFunction, constructible_from_function_rvalueref)
{
  const GenericFunction<d, r> some_func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  GridFunction<E, r> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(VectorGridFunction, convertible_from_function_rvalueref)
{
  const GenericFunction<d, r> func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  auto nm = accepts_vector_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(VectorGridFunction, constructible_from_function_lvalueptr)
{
  GridFunction<E, r> func{new GenericFunction<d, r>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(VectorGridFunction, convertible_from_function_lvalueptr)
{
  accepts_vector_grid_function(new GenericFunction<d, r>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"});
}
GTEST_TEST(VectorGridFunction, constructible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r> some_func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  GridFunction<E, r> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(VectorGridFunction, convertible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r> func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  auto nm = accepts_vector_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(VectorGridFunction, constructible_from_grid_function_lvalueptr)
{
  GridFunction<E, r> func{
      new GenericGridFunction<E, r>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(VectorGridFunction, convertible_from_grid_function_lvalueptr)
{
  accepts_vector_grid_function(
      new GenericGridFunction<E, r>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"});
}


// SquareMatrixGridFunction

GTEST_TEST(SquareMatrixGridFunction, constructible_from_double)
{
  GridFunction<E, r, r> DXTC_UNUSED(func){1.};
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_double)
{
  accepts_square_matrix_grid_function(1.);
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_FieldMatrix)
{
  GridFunction<E, r, r> DXTC_UNUSED(func){FieldMatrix<double, r, r>()};
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_FieldMatrix)
{
  accepts_square_matrix_grid_function(FieldMatrix<double, r, r>());
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_XT_Common_FieldMatrix)
{
  GridFunction<E, r, r> DXTC_UNUSED(func){XT::Common::FieldMatrix<double, r, r>()};
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_XT_Common_FieldMatrix)
{
  accepts_square_matrix_grid_function(XT::Common::FieldMatrix<double, r, r>());
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_function_rvalueref)
{
  const GenericFunction<d, r, r> some_func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  GridFunction<E, r, r> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_function_rvalueref)
{
  const GenericFunction<d, r, r> func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  auto nm = accepts_square_matrix_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_function_lvalueptr)
{
  GridFunction<E, r, r> func{new GenericFunction<d, r, r>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_function_lvalueptr)
{
  accepts_square_matrix_grid_function(new GenericFunction<d, r, r>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"});
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r, r> some_func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  GridFunction<E, r, r> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r, r> func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  auto nm = accepts_square_matrix_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(SquareMatrixGridFunction, constructible_from_grid_function_lvalueptr)
{
  GridFunction<E, r, r> func{
      new GenericGridFunction<E, r, r>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(SquareMatrixGridFunction, convertible_from_grid_function_lvalueptr)
{
  accepts_square_matrix_grid_function(
      new GenericGridFunction<E, r, r>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"});
}


// MatrixGridFunction

GTEST_TEST(MatrixGridFunction, constructible_from_FieldMatrix)
{
  GridFunction<E, r, rC> DXTC_UNUSED(func){FieldMatrix<double, r, rC>()};
}
GTEST_TEST(MatrixGridFunction, convertible_from_FieldMatrix)
{
  accepts_matrix_grid_function(FieldMatrix<double, r, rC>());
}
GTEST_TEST(MatrixGridFunction, constructible_from_XT_Common_FieldMatrix)
{
  GridFunction<E, r, rC> DXTC_UNUSED(func){XT::Common::FieldMatrix<double, r, rC>()};
}
GTEST_TEST(MatrixGridFunction, convertible_from_XT_Common_FieldMatrix)
{
  accepts_matrix_grid_function(XT::Common::FieldMatrix<double, r, rC>());
}
GTEST_TEST(MatrixGridFunction, constructible_from_function_rvalueref)
{
  const GenericFunction<d, r, rC> some_func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  GridFunction<E, r, rC> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(MatrixGridFunction, convertible_from_function_rvalueref)
{
  const GenericFunction<d, r, rC> func{0, [](auto&, auto&) { return 1.; }, "THE_NAME"};
  auto nm = accepts_matrix_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(MatrixGridFunction, constructible_from_function_lvalueptr)
{
  GridFunction<E, r, rC> func{new GenericFunction<d, r, rC>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(MatrixGridFunction, convertible_from_function_lvalueptr)
{
  accepts_matrix_grid_function(new GenericFunction<d, r, rC>{0, [](auto&, auto&) { return 1.; }, "THE_NAME"});
}
GTEST_TEST(MatrixGridFunction, constructible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r, rC> some_func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  GridFunction<E, r, rC> func{some_func};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(MatrixGridFunction, convertible_from_grid_function_rvalueref)
{
  const GenericGridFunction<E, r, rC> func{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"};
  auto nm = accepts_matrix_grid_function(func);
  EXPECT_EQ(std::string("THE_NAME"), nm);
}
GTEST_TEST(MatrixGridFunction, constructible_from_grid_function_lvalueptr)
{
  GridFunction<E, r, rC> func{
      new GenericGridFunction<E, r, rC>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"}};
  EXPECT_EQ(std::string("THE_NAME"), func.name());
}
GTEST_TEST(MatrixGridFunction, convertible_from_grid_function_lvalueptr)
{
  accepts_matrix_grid_function(
      new GenericGridFunction<E, r, rC>{0, [](auto&) {}, [](auto&, auto&) { return 1.; }, {}, "THE_NAME"});
}
