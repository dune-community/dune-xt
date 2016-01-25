// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2014, 2016)
//   Rene Milk       (2014 - 2015)
//   Tobias Leibner  (2014 - 2015)

#include <dune/xt/common/test/main.hxx>

#include "container.hh"

using namespace Dune;
using namespace Dune::XT;

static const size_t dim = 4;

struct ContainerTest : public ::testing::Test
{
  void fulfills_interface() const
  {
    typedef TESTCONTAINERTYPE ContainerImp;
    // static tests
    typedef typename ContainerImp::Traits Traits;
    // * of the traits
    typedef typename Traits::derived_type T_derived_type;
    static_assert(std::is_same<ContainerImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    typedef typename Traits::ScalarType T_ScalarType;
    typedef typename Traits::RealType T_RealType;
    // * of the container as itself (aka the derived type)
    typedef typename ContainerImp::ScalarType D_ScalarType;
    typedef typename ContainerImp::RealType D_RealType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    static_assert(std::is_same<T_RealType, D_RealType>::value, "RealType of derived_type has to be the correct Type!");
    // * of the container as the interface
    typedef typename XT::LA::ContainerInterface<Traits, D_ScalarType> InterfaceType;
    typedef typename InterfaceType::derived_type I_derived_type;
    typedef typename InterfaceType::ScalarType I_ScalarType;
    typedef typename InterfaceType::RealType I_RealType;
    static_assert(std::is_same<ContainerImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    static_assert(std::is_same<T_RealType, I_RealType>::value,
                  "ScalarType of derived_type has to be the correct Type!");

    // dynamic tests
    // * of the container as itself (aka the derived type)
    ContainerImp DUNE_UNUSED(d_empty);
    ContainerImp d_by_size = ContainerFactory<ContainerImp>::create(dim);
    ContainerImp d_copy_constructor(d_by_size);
    ContainerImp DUNE_UNUSED(d_copy_assignment) = d_by_size;
    ContainerImp d_deep_copy = d_by_size.copy();
    d_by_size.scal(D_ScalarType(1));
    d_by_size.axpy(D_ScalarType(1), d_deep_copy);
    EXPECT_TRUE(d_by_size.has_equal_shape(d_deep_copy));
    // * of the container as the interface
    InterfaceType& i_by_size = static_cast<InterfaceType&>(d_by_size);
    ContainerImp i_deep_copy = i_by_size.copy();
    i_by_size.scal(I_ScalarType(1));
    i_by_size.axpy(I_ScalarType(1), i_deep_copy);
  } // void fulfills_interface() const
}; // struct ContainerTest

TEST_F(ContainerTest, fulfills_interface)
{
  this->fulfills_interface();
}
