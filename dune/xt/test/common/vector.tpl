// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2020)

#include <dune/xt/test/main.hxx>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/vector.hh>

#include <dune/xt/test/la/container.hh>

using namespace Dune;
using namespace Dune::XT;
using namespace Dune::XT::Common;

static constexpr size_t dim = 4;

{% for T_NAME, V_TYPE in config.testtypes %}
struct VectorTest_{{T_NAME}} : public ::testing::Test
{
  using VectorImp = {{V_TYPE}};
  using V = VectorAbstraction<VectorImp>;
  using ScalarType = typename V::ScalarType;

  void test_free_functions() const
  {
    auto vec = ContainerFactory<VectorImp>::create(dim);
    for (size_t ii = 0; ii < dim; ++ii)
      vec[ii] = static_cast<ScalarType>(ii);

    // test sup_norm
    EXPECT_EQ(sup_norm(vec), static_cast<ScalarType>(3));

    // test zeros_like
    auto zeros = zeros_like(vec);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(zeros, XT::Common::FieldVector<ScalarType, dim>(0)));

    // test create
    auto ones = create<VectorImp>(dim, static_cast<ScalarType>(1));
    EXPECT_TRUE(XT::Common::FloatCmp::eq(ones, XT::Common::FieldVector<ScalarType, dim>(1)));

    // test data
    [[maybe_unused]] ScalarType* data_ptr = XT::Common::data(vec);
    const VectorImp const_vec = vec;
    [[maybe_unused]] const ScalarType* const_data_ptr = XT::Common::data(const_vec);
  } // void test_free_functions() const
}; // struct VectorTest

TEST_F(VectorTest_{{T_NAME}}, test_free_functions)
{
  this->test_free_functions();
}
{% endfor %}
