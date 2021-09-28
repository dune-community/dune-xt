// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2018)
//   René Fritze     (2014 - 2017, 2019)
//   Tobias Leibner  (2014 - 2017, 2020)

#include <dune/xt/test/main.hxx>

#include <dune/xt/common/float_cmp.hh>

#include <dune/xt/test/la/container.hh>

using namespace Dune;
using namespace Dune::XT;

static constexpr size_t dim = 4;
{% for NAME,CIMP in config.container %}
struct ContainerTest{{NAME}} : public ::testing::Test
{
  void fulfills_interface() const
  {
    // static tests
    using ContainerImp = {{CIMP}};
    using Traits = typename ContainerImp::Traits;
    // * of the traits
    using T_derived_type = typename Traits::derived_type;
    static_assert(std::is_same<ContainerImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    using T_ScalarType = typename Traits::ScalarType;
    // * of the container as itself (aka the derived type)
    using D_ScalarType = typename ContainerImp::ScalarType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // * of the container as the interface
    using InterfaceType = typename XT::LA::ContainerInterface<Traits, D_ScalarType>;
    using I_derived_type = typename InterfaceType::derived_type;
    using I_ScalarType = typename InterfaceType::ScalarType;
    static_assert(std::is_same<ContainerImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");

    // dynamic tests
    // * of the container as itself (aka the derived type)
    [[maybe_unused]] ContainerImp d_empty;
    ContainerImp d_by_size = ContainerFactory<ContainerImp>::create(dim);
    ContainerImp d_copy_constructor(d_by_size);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(d_by_size, d_copy_constructor)); // checks that copy works
    d_copy_constructor.scal(D_ScalarType(10));
    EXPECT_TRUE(XT::Common::FloatCmp::ne(d_by_size, d_copy_constructor)); // checks that only copied container changes
    ContainerImp d_copy_assignment = d_by_size;
    EXPECT_TRUE(XT::Common::FloatCmp::eq(d_by_size, d_copy_assignment));
    d_copy_assignment.scal(D_ScalarType(10));
    EXPECT_TRUE(XT::Common::FloatCmp::ne(d_by_size, d_copy_assignment));
    ContainerImp d_deep_copy = d_by_size.copy();
    EXPECT_TRUE(XT::Common::FloatCmp::eq(d_by_size, d_deep_copy));
    d_deep_copy.scal(D_ScalarType(10));
    EXPECT_TRUE(XT::Common::FloatCmp::ne(d_by_size, d_deep_copy));
    d_by_size.scal(D_ScalarType(1));
    d_by_size.axpy(D_ScalarType(1), d_deep_copy);
    EXPECT_TRUE(d_by_size.has_equal_shape(d_deep_copy));
    // * of the container as the interface
    InterfaceType& i_by_size = static_cast<InterfaceType&>(d_by_size);
    ContainerImp i_deep_copy = i_by_size.copy();
    i_by_size.scal(I_ScalarType(1));
    i_by_size.axpy(I_ScalarType(1), i_deep_copy);
    i_by_size = d_copy_constructor;
    EXPECT_TRUE(XT::Common::FloatCmp::eq(d_by_size, d_copy_constructor));
  } // void fulfills_interface() const
}; // struct ContainerTest

TEST_F(ContainerTest{{NAME}}, fulfills_interface)
{
  this->fulfills_interface();
}
{% endfor %}
