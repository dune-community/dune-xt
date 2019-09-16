// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2017)

#include <dune/xt/common/test/main.hxx>
#include <dune/xt/common/vector.hh>

#include <dune/xt/la/test/container.hh>

using namespace Dune;
using namespace Dune::XT;

static const size_t dim = 4;

{% for T_NAME, V_TYPE in config.testtypes %}
struct VectorTest_{{T_NAME}} : public ::testing::Test
{
  typedef {{V_TYPE}} VectorImp;

  void fulfills_interface() const
  {
    // static tests
    typedef typename VectorImp::Traits Traits;
    // * of the traits
    typedef typename Traits::derived_type T_derived_type;
    static_assert(std::is_same<VectorImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    typedef typename Traits::ScalarType T_ScalarType;
    typedef typename Traits::RealType T_RealType;
    // * of the vector as itself (aka the derived type)
    typedef typename VectorImp::ScalarType D_ScalarType;
    typedef typename VectorImp::RealType D_RealType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    static_assert(std::is_same<T_RealType, D_RealType>::value, "RealType of derived_type has to be the correct Type!");
    // * of the vector as the interface
    typedef typename XT::LA::VectorInterface<Traits, D_ScalarType> InterfaceType;
    typedef typename InterfaceType::derived_type I_derived_type;
    typedef typename InterfaceType::ScalarType I_ScalarType;
    typedef typename InterfaceType::RealType I_RealType;
    static_assert(std::is_same<VectorImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    static_assert(std::is_same<T_RealType, I_RealType>::value, "RealType of derived_type has to be the correct Type!");
    // dynamic tests
    // * of the vector as itself (aka the derived type)
    VectorImp d_by_size(dim);
    VectorImp d_by_size_and_value(dim, D_ScalarType(0));
    size_t d_size = d_by_size.size();
    EXPECT_EQ(dim, d_size);
    for (size_t ii = 0; ii < d_size; ++ii) {
      d_by_size_and_value.set_entry(ii, D_ScalarType(0.5) + D_ScalarType(ii));
      d_by_size_and_value.add_to_entry(ii, D_ScalarType(0.5) + D_ScalarType(ii));
      EXPECT_FALSE(
          Common::FloatCmp::ne(d_by_size_and_value.get_entry(ii), D_ScalarType(2) * D_ScalarType(ii) + D_ScalarType(1)))
          << d_by_size_and_value.get_entry(ii);
      EXPECT_FALSE(Common::FloatCmp::ne(d_by_size_and_value.get_entry(ii), d_by_size_and_value.get_entry(ii)))
          << d_by_size_and_value.get_entry(ii);
    }
    EXPECT_TRUE(d_by_size.almost_equal(d_by_size));
    d_by_size_and_value.scal(D_ScalarType(0));
    EXPECT_TRUE(d_by_size_and_value.almost_equal(d_by_size));
    D_ScalarType d_dot = d_by_size.dot(d_by_size_and_value);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(D_RealType(0), d_dot);
    D_ScalarType d_l1_norm = d_by_size.l1_norm();
    D_ScalarType d_l1_norm2 = d_by_size.one_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(D_RealType(0), d_l1_norm);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(d_l1_norm, d_l1_norm2);
    D_ScalarType d_l2_norm = d_by_size.l2_norm();
    D_ScalarType d_l2_norm2 = d_by_size.two_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(D_RealType(0), d_l2_norm);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(d_l2_norm, d_l2_norm2);
    D_ScalarType d_sup_norm = d_by_size.sup_norm();
    D_ScalarType d_sup_norm2 = d_by_size.inf_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(D_RealType(0), d_sup_norm);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(d_sup_norm, d_sup_norm2);
    VectorImp d_ones(dim, D_ScalarType(1));
    std::pair<size_t, D_ScalarType> d_amax = d_ones.amax();
    EXPECT_EQ(0, d_amax.first);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(D_RealType(1), d_amax.second);
    d_ones.add(d_by_size, d_by_size_and_value);
    EXPECT_TRUE(d_by_size_and_value.almost_equal(d_ones));
    VectorImp d_added = d_ones.add(d_by_size);
    EXPECT_TRUE(d_added.almost_equal(d_by_size_and_value));
    d_added.iadd(d_by_size);
    EXPECT_TRUE(d_added.almost_equal(d_ones));
    d_ones.sub(d_by_size, d_by_size_and_value);
    EXPECT_TRUE(d_by_size_and_value.almost_equal(d_ones));
    VectorImp d_subtracted = d_ones.sub(d_by_size);
    EXPECT_TRUE(d_subtracted.almost_equal(d_by_size_and_value));
    d_subtracted.isub(d_by_size);
    EXPECT_TRUE(d_subtracted.almost_equal(d_ones));
    // * of the vector as the interface
    VectorImp d_by_size_2(dim);
    VectorImp d_by_size_and_value_2(dim, D_ScalarType(1));
    InterfaceType& i_by_size = static_cast<InterfaceType&>(d_by_size_2);
    InterfaceType& i_by_size_and_value = static_cast<InterfaceType&>(d_by_size_and_value_2);
    EXPECT_TRUE(i_by_size.almost_equal(d_by_size_2));
    i_by_size_and_value.scal(I_ScalarType(0));
    EXPECT_TRUE(i_by_size_and_value.almost_equal(d_by_size_2));
    I_ScalarType i_dot = i_by_size.dot(d_by_size_and_value_2);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(I_RealType(0), i_dot);
    I_ScalarType i_l1_norm = i_by_size.l1_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(I_RealType(0), i_l1_norm);
    I_ScalarType i_l2_norm = i_by_size.l2_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(I_RealType(0), i_l2_norm);
    I_ScalarType i_sup_norm = i_by_size.sup_norm();
    EXPECT_DOUBLE_OR_COMPLEX_EQ(I_RealType(0), i_sup_norm);
    VectorImp i_ones(dim, I_ScalarType(1));
    std::pair<size_t, I_ScalarType> i_amax = i_ones.amax();
    EXPECT_EQ(0, i_amax.first);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(I_RealType(1), i_amax.second);
    i_ones.add(d_by_size_2, d_by_size_and_value_2);
    EXPECT_TRUE(i_by_size_and_value.almost_equal(i_ones));
    VectorImp i_added = i_ones.add(d_by_size_2);
    EXPECT_TRUE(i_added.almost_equal(d_by_size_and_value_2));
    i_added.iadd(d_by_size_2);
    EXPECT_TRUE(i_added.almost_equal(i_ones));
    i_ones.sub(d_by_size_2, d_by_size_and_value_2);
    EXPECT_TRUE(i_by_size_and_value.almost_equal(i_ones));
    VectorImp i_subtracted = i_ones.sub(d_by_size_2);
    EXPECT_TRUE(i_subtracted.almost_equal(d_by_size_and_value_2));
    i_subtracted.isub(d_by_size_2);
    EXPECT_TRUE(i_subtracted.almost_equal(i_ones));
  } // void fulfills_interface() const

  void produces_correct_results() const
  {
    typedef typename VectorImp::ScalarType ScalarType;
    typedef typename VectorImp::RealType RealType;
    using V = XT::Common::VectorAbstraction<VectorImp>;

    // create test vectors
    VectorImp zeros(dim); // [0, 0, 0, 0]
    VectorImp ones(dim, ScalarType(1)); // [1, 1, 1, 1]
    VectorImp countingup(dim); //[0, 1, 2, 3]
    for (size_t ii = 0; ii < dim; ++ii)
      countingup.set_entry(ii, ScalarType(ii));
    VectorImp testvector_1(dim); // [0, -2, 2, 1]
    testvector_1.set_entry(0, ScalarType(0));
    testvector_1.set_entry(1, ScalarType(-2));
    testvector_1.set_entry(2, ScalarType(2));
    testvector_1.set_entry(3, ScalarType(1));
    VectorImp testvector_2(dim); // [0, 2, -2, 1]
    testvector_2.set_entry(0, ScalarType(0));
    testvector_2.set_entry(1, ScalarType(2));
    testvector_2.set_entry(2, ScalarType(-2));
    testvector_2.set_entry(3, ScalarType(1));
    VectorImp testvector_3(dim); // [-1, 1, -1, 1]
    testvector_3.set_entry(0, ScalarType(-1));
    testvector_3.set_entry(1, ScalarType(1));
    testvector_3.set_entry(2, ScalarType(-1));
    testvector_3.set_entry(3, ScalarType(1));
    VectorImp testvector_4(dim); // [0, 3, -2, 0]
    testvector_4.set_entry(0, ScalarType(0));
    testvector_4.set_entry(1, ScalarType(3));
    testvector_4.set_entry(2, ScalarType(-2));
    testvector_4.set_entry(3, ScalarType(0));
    VectorImp testvector_5(dim); // [1.25, 0, 2.5, -3.5]
    testvector_5.set_entry(0, ScalarType(1.25));
    testvector_5.set_entry(1, ScalarType(0));
    testvector_5.set_entry(2, ScalarType(2.5));
    testvector_5.set_entry(3, ScalarType(-3.5));

    // test resize, should not change existing entries
    if (!V::has_static_size) {
        try {
        auto testvector_1_saved = testvector_1;
        testvector_1.resize(dim-1);
        EXPECT_EQ(dim-1, testvector_1.size());
        EXPECT_EQ(dim, testvector_1_saved.size());
        for (size_t ii = 0; ii < dim-1; ++ii)
          EXPECT_DOUBLE_OR_COMPLEX_EQ(std::real(testvector_1_saved.get_entry(ii)), testvector_1.get_entry(ii));
        testvector_1.resize(dim);
        testvector_1.set_entry(dim-1, testvector_1_saved[dim-1]);
        for (size_t ii = 0; ii < dim; ++ii)
          EXPECT_DOUBLE_OR_COMPLEX_EQ(std::real(testvector_1_saved.get_entry(ii)), testvector_1.get_entry(ii));
        } catch(XT::Common::Exceptions::you_are_using_this_wrong) {
        }
    }

    // test get_entry()
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), testvector_1.get_entry(0));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(-2), testvector_1.get_entry(1));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(2), testvector_1.get_entry(2));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1), testvector_1.get_entry(3));

    // test operator[]
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), testvector_1[0]);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(-2), testvector_1[1]);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(2), testvector_1[2]);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1), testvector_1[3]);
    VectorImp testvector_1_copy = testvector_1;
    const ScalarType& entry0 = testvector_1_copy[0];
    ScalarType& entry3 = testvector_1_copy[3];
    testvector_1.scal(ScalarType(2));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), entry0);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1), entry3);
    testvector_1.scal(ScalarType(0.5));
    testvector_1_copy.scal(ScalarType(3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(0), entry0);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(3), entry3);
    entry3 = ScalarType(42);
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(42), testvector_1_copy.get_entry(3));
    EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(1), testvector_1.get_entry(3));

    // test iterators
    const auto it_end = countingup.end();
    size_t count = 0;
    for (auto it = countingup.begin(); it != it_end; ++it)
      EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(count++), *it);
    count = 0;
    for (const auto& entry : countingup)
      EXPECT_DOUBLE_OR_COMPLEX_EQ(RealType(count++), entry);

    // test amax()
    std::pair<size_t, RealType> amax = zeros.amax();
    EXPECT_EQ(0, amax.first);
    EXPECT_DOUBLE_EQ(RealType(0), amax.second);
    amax = ones.amax();
    EXPECT_EQ(0, amax.first);
    EXPECT_DOUBLE_EQ(RealType(1), amax.second);
    amax = countingup.amax();
    EXPECT_EQ(3, amax.first);
    EXPECT_DOUBLE_EQ(RealType(3), amax.second);
    amax = testvector_1.amax();
    EXPECT_EQ(1, amax.first);
    EXPECT_DOUBLE_EQ(RealType(2), amax.second);
    amax = testvector_2.amax();
    EXPECT_EQ(1, amax.first);
    EXPECT_DOUBLE_EQ(RealType(2), amax.second);
    amax = testvector_3.amax();
    EXPECT_EQ(0, amax.first);
    EXPECT_DOUBLE_EQ(RealType(1), amax.second);
    amax = testvector_4.amax();
    EXPECT_EQ(1, amax.first);
    EXPECT_DOUBLE_EQ(RealType(3), amax.second);
    amax = testvector_5.amax();
    EXPECT_EQ(3, amax.first);
    EXPECT_DOUBLE_EQ(RealType(3.5), amax.second);

    // test max, min
    if (!std::is_arithmetic<ScalarType>::value) {
      EXPECT_THROW(ones.min(), XT::LA::Exceptions::not_available);
      EXPECT_THROW(ones.max(), XT::LA::Exceptions::not_available);
    } else {
      EXPECT_EQ(-2., testvector_4.min());
      EXPECT_EQ(3., testvector_4.max());
      EXPECT_EQ(-3.5, testvector_5.min());
      EXPECT_EQ(2.5, testvector_5.max());
    }

    // test l1_norm()
    RealType l1_norm = zeros.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(0), l1_norm);
    l1_norm = ones.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(4), l1_norm);
    l1_norm = countingup.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(6), l1_norm);
    l1_norm = testvector_1.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(5), l1_norm);
    l1_norm = testvector_2.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(5), l1_norm);
    l1_norm = testvector_3.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(4), l1_norm);
    l1_norm = testvector_4.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(5), l1_norm);
    l1_norm = testvector_5.l1_norm();
    EXPECT_DOUBLE_EQ(RealType(7.25), l1_norm);

    // test l2_norm()
    RealType l2_norm = zeros.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(0), l2_norm);
    l2_norm = ones.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(2), l2_norm);
    l2_norm = countingup.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(std::sqrt(14)), l2_norm);
    l2_norm = testvector_1.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(3), l2_norm);
    l2_norm = testvector_2.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(3), l2_norm);
    l2_norm = testvector_3.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(2), l2_norm);
    l2_norm = testvector_4.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(std::sqrt(13)), l2_norm);
    l2_norm = testvector_5.l2_norm();
    EXPECT_DOUBLE_EQ(RealType(std::sqrt(20.0625)), l2_norm);

    // test sup_norm()
    RealType sup_norm = zeros.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(0), sup_norm);
    sup_norm = ones.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(1), sup_norm);
    sup_norm = countingup.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(3), sup_norm);
    sup_norm = testvector_1.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(2), sup_norm);
    sup_norm = testvector_2.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(2), sup_norm);
    sup_norm = testvector_3.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(1), sup_norm);
    sup_norm = testvector_4.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(3), sup_norm);
    sup_norm = testvector_5.sup_norm();
    EXPECT_DOUBLE_EQ(RealType(3.5), sup_norm);

    // test dot(), operator*
    ScalarType dot = ones.dot(zeros);
    ScalarType dot_operator = ones * zeros;
    ScalarType dot2 = zeros.dot(ones);
    ScalarType dot_operator_2 = zeros * ones;
    EXPECT_TRUE(Common::FloatCmp::eq(dot, ScalarType(0)) && Common::FloatCmp::eq(dot, dot2)
                && Common::FloatCmp::eq(dot_operator, dot_operator_2)
                && Common::FloatCmp::eq(dot, dot_operator))
        << "These should all equal 0: " << dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2;
    dot = ones.dot(ones);
    dot_operator = ones * ones;
    EXPECT_TRUE(Common::FloatCmp::eq(dot, ScalarType(4)) && Common::FloatCmp::eq(dot_operator, ScalarType(4)))
        << "These should equal 4: " << dot << ", " << dot_operator;
    dot = ones.dot(testvector_3);
    dot_operator = ones * testvector_3;
    dot2 = testvector_3.dot(ones);
    dot_operator_2 = testvector_3 * ones;
    EXPECT_TRUE(Common::FloatCmp::eq(dot, ScalarType(0)) && Common::FloatCmp::eq(dot, dot2)
                && Common::FloatCmp::eq(dot_operator, dot_operator_2)
                && Common::FloatCmp::eq(dot, dot_operator))
        << "These should all equal 0: " << dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2;
    dot = countingup.dot(testvector_5);
    dot_operator = countingup * testvector_5;
    dot2 = testvector_5.dot(countingup);
    dot_operator_2 = testvector_5 * countingup;
    EXPECT_TRUE(Common::FloatCmp::eq(dot, ScalarType(-5.5)) && Common::FloatCmp::eq(dot, dot2)
                && Common::FloatCmp::eq(dot_operator, dot_operator_2)
                && Common::FloatCmp::eq(dot, dot_operator))
        << "These should all equal -5.5: " << dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2;
    dot = testvector_3.dot(testvector_5);
    dot_operator = testvector_3 * testvector_5;
    dot2 = testvector_5.dot(testvector_3);
    dot_operator_2 = testvector_5 * testvector_3;
    EXPECT_TRUE(Common::FloatCmp::eq(dot, ScalarType(-7.25)) && Common::FloatCmp::eq(dot, dot2)
                && Common::FloatCmp::eq(dot_operator, dot_operator_2)
                && Common::FloatCmp::eq(dot, dot_operator))
        << "These should all equal -7.25: " << dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2;

    // test operator==
    EXPECT_TRUE(zeros == zeros);
    EXPECT_TRUE(ones == ones);
    EXPECT_TRUE(countingup == countingup);
    EXPECT_TRUE(testvector_1 == testvector_1);
    EXPECT_TRUE(testvector_2 == testvector_2);
    EXPECT_TRUE(testvector_3 == testvector_3);
    EXPECT_TRUE(testvector_4 == testvector_4);
    EXPECT_TRUE(testvector_5 == testvector_5);

    // test operator!=
    EXPECT_TRUE(zeros != ones);
    EXPECT_TRUE(testvector_1 != testvector_2);
    EXPECT_TRUE(testvector_3 != testvector_4);
    EXPECT_TRUE(testvector_3 != ones);
    EXPECT_TRUE(testvector_4 != testvector_5);
    EXPECT_FALSE(ones != ones);
    EXPECT_FALSE(testvector_3 != testvector_3);

    // test scal, operator*=
    VectorImp scaled = zeros;
    VectorImp scaled_by_operator = zeros;
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    EXPECT_EQ(zeros, scaled);
    EXPECT_EQ(zeros, scaled_by_operator);
    scaled = ones;
    scaled_by_operator = ones;
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    EXPECT_EQ(ones, scaled);
    EXPECT_EQ(ones, scaled_by_operator);
    scaled.scal(ScalarType(0));
    scaled_by_operator *= ScalarType(0);
    EXPECT_EQ(zeros, scaled);
    EXPECT_EQ(zeros, scaled_by_operator);
    scaled = testvector_1;
    scaled.scal(ScalarType(2));
    scaled_by_operator = testvector_1;
    scaled_by_operator *= ScalarType(2);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(2) * testvector_1.get_entry(ii), scaled.get_entry(ii)));
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(2) * testvector_1.get_entry(ii), scaled_by_operator.get_entry(ii)));
    }
    scaled = testvector_3;
    scaled.scal(ScalarType(-2));
    scaled_by_operator = testvector_3;
    scaled_by_operator *= ScalarType(-2);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(-2) * testvector_3.get_entry(ii), scaled.get_entry(ii)));
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(-2) * testvector_3.get_entry(ii), scaled_by_operator.get_entry(ii)));
    }
    scaled = countingup;
    scaled.scal(ScalarType(2.2));
    scaled_by_operator = countingup;
    scaled_by_operator *= ScalarType(2.2);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(2.2) * countingup.get_entry(ii), scaled.get_entry(ii)));
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(2.2) * countingup.get_entry(ii), scaled_by_operator.get_entry(ii)));
    }
    scaled = testvector_5;
    scaled.scal(ScalarType(-3.75));
    scaled_by_operator = testvector_5;
    scaled_by_operator *= ScalarType(-3.75);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(-3.75) * testvector_5.get_entry(ii), scaled.get_entry(ii)));
      EXPECT_TRUE(
          Common::FloatCmp::eq(ScalarType(-3.75) * testvector_5.get_entry(ii), scaled_by_operator.get_entry(ii)));
    }
    VectorImp a = ones;
    a.scal(ScalarType(0));
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }

    // test operator+, operator+=, add, iadd
    VectorImp sum_operator_plus = zeros + ones;
    VectorImp sum_operator_iplus = zeros;
    sum_operator_iplus += ones;
    VectorImp sum_add_1 = zeros.add(ones);
    VectorImp sum_add_2(dim);
    zeros.add(ones, sum_add_2);
    VectorImp sum_iadd = zeros;
    sum_iadd.iadd(ones);
    VectorImp sum_correct(dim);
    sum_correct.set_entry(0, ScalarType(1));
    sum_correct.set_entry(1, ScalarType(1));
    sum_correct.set_entry(2, ScalarType(1));
    sum_correct.set_entry(3, ScalarType(1));
    EXPECT_FALSE(sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
                 || sum_add_1 != sum_operator_plus
                 || sum_add_1 != sum_correct)
        << "These should all be equal:\n"
        << sum_add_1 << ",\n"
        << sum_add_2 << ",\n"
        << sum_operator_iplus << ",\n"
        << sum_operator_plus << ",\n"
        << sum_iadd << ",\n"
        << sum_correct;

    sum_operator_plus = countingup + testvector_1;
    sum_operator_iplus = countingup;
    sum_operator_iplus += testvector_1;
    sum_add_1 = countingup.add(testvector_1);
    countingup.add(testvector_1, sum_add_2);
    sum_iadd = countingup;
    sum_iadd.iadd(testvector_1);
    sum_correct.set_entry(0, ScalarType(0));
    sum_correct.set_entry(1, ScalarType(-1));
    sum_correct.set_entry(2, ScalarType(4));
    sum_correct.set_entry(3, ScalarType(4));
    EXPECT_FALSE(sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
                 || sum_add_1 != sum_operator_plus
                 || sum_add_1 != sum_correct)
        << "These should all be equal:\n"
        << sum_add_1 << ",\n"
        << sum_add_2 << ",\n"
        << sum_operator_iplus << ",\n"
        << sum_operator_plus << ",\n"
        << sum_iadd << ",\n"
        << sum_correct;

    sum_operator_plus = testvector_3 + testvector_5;
    sum_operator_iplus = testvector_3;
    sum_operator_iplus += testvector_5;
    sum_add_1 = testvector_3.add(testvector_5);
    testvector_3.add(testvector_5, sum_add_2);
    sum_iadd = testvector_3;
    sum_iadd.iadd(testvector_5);
    sum_correct.set_entry(0, ScalarType(0.25));
    sum_correct.set_entry(1, ScalarType(1));
    sum_correct.set_entry(2, ScalarType(1.5));
    sum_correct.set_entry(3, ScalarType(-2.5));
    EXPECT_FALSE(sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
                 || sum_add_1 != sum_operator_plus
                 || sum_add_1 != sum_correct)
        << "These should all be equal:\n"
        << sum_add_1 << ",\n"
        << sum_add_2 << ",\n"
        << sum_operator_iplus << ",\n"
        << sum_operator_plus << ",\n"
        << sum_iadd << ",\n"
        << sum_correct;

    a = ones;
    a += testvector_3;
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }
    a = ones;
    a.iadd(testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }
    a = ones;
    ones.add(testvector_3, a);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }

    // test operator-, operator-=, sub, isub
    VectorImp diff_operator_minus = zeros - ones;
    VectorImp diff_operator_iminus = zeros;
    diff_operator_iminus -= ones;
    VectorImp diff_sub_1 = zeros.sub(ones);
    VectorImp diff_sub_2(dim);
    zeros.sub(ones, diff_sub_2);
    VectorImp diff_isub = zeros;
    diff_isub.isub(ones);
    VectorImp diff_correct(dim);
    diff_correct.set_entry(0, ScalarType(-1));
    diff_correct.set_entry(1, ScalarType(-1));
    diff_correct.set_entry(2, ScalarType(-1));
    diff_correct.set_entry(3, ScalarType(-1));
    EXPECT_TRUE(diff_operator_iminus == diff_operator_minus && diff_sub_1 == diff_sub_2 && diff_isub == diff_sub_1
                && diff_sub_1 == diff_operator_minus
                && diff_sub_1 == diff_correct)
        << "These should all be equal:\n"
        << diff_sub_1 << ",\n"
        << diff_sub_2 << ",\n"
        << diff_operator_iminus << ",\n"
        << diff_operator_minus << ",\n"
        << diff_isub << ",\n"
        << diff_correct;

    diff_operator_minus = testvector_1 - testvector_4;
    diff_operator_iminus = testvector_1;
    diff_operator_iminus -= testvector_4;
    diff_sub_1 = testvector_1.sub(testvector_4);
    testvector_1.sub(testvector_4, diff_sub_2);
    diff_isub = testvector_1;
    diff_isub.isub(testvector_4);
    diff_correct.set_entry(0, ScalarType(-0));
    diff_correct.set_entry(1, ScalarType(-5));
    diff_correct.set_entry(2, ScalarType(4));
    diff_correct.set_entry(3, ScalarType(1));
    EXPECT_TRUE(diff_operator_iminus == diff_operator_minus && diff_sub_1 == diff_sub_2 && diff_isub == diff_sub_1
                && diff_sub_1 == diff_operator_minus
                && diff_sub_1 == diff_correct)
        << "These should all be equal:\n"
        << diff_sub_1 << ",\n"
        << diff_sub_2 << ",\n"
        << diff_operator_iminus << ",\n"
        << diff_operator_minus << ",\n"
        << diff_isub << ",\n"
        << diff_correct;

    diff_operator_minus = testvector_5 - testvector_2;
    diff_operator_iminus = testvector_5;
    diff_operator_iminus -= testvector_2;
    diff_sub_1 = testvector_5.sub(testvector_2);
    testvector_5.sub(testvector_2, diff_sub_2);
    diff_isub = testvector_5;
    diff_isub.isub(testvector_2);
    diff_correct.set_entry(0, ScalarType(1.25));
    diff_correct.set_entry(1, ScalarType(-2));
    diff_correct.set_entry(2, ScalarType(4.5));
    diff_correct.set_entry(3, ScalarType(-4.5));
    EXPECT_TRUE(diff_operator_iminus == diff_operator_minus && diff_sub_1 == diff_sub_2 && diff_isub == diff_sub_1
                && diff_sub_1 == diff_operator_minus
                && diff_sub_1 == diff_correct)
        << "These should all be equal:\n"
        << diff_sub_1 << ",\n"
        << diff_sub_2 << ",\n"
        << diff_operator_iminus << ",\n"
        << diff_operator_minus << ",\n"
        << diff_isub << ",\n"
        << diff_correct;

    a = ones;
    a -= testvector_3;
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }
    a = ones;
    a.isub(testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }
    a = ones;
    ones.sub(testvector_3, a);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }

    // test operator= for scalars
    VectorImp result_eq = countingup;
    result_eq = ScalarType(0);
    EXPECT_EQ(zeros, result_eq);
    result_eq = ScalarType(-2.25);
    VectorImp result_eq_expected = ones;
    result_eq_expected *= ScalarType(-2.25);
    EXPECT_EQ(result_eq_expected, result_eq);

    // test axpy
    VectorImp result_axpy = zeros;
    result_axpy.axpy(ScalarType(1), ones);
    EXPECT_EQ(ones, result_axpy);
    result_axpy = testvector_5;
    result_axpy.axpy(ScalarType(2.75), testvector_3);
    VectorImp correct_result = testvector_3;
    correct_result.scal(ScalarType(2.75));
    correct_result += testvector_5;
    EXPECT_EQ(correct_result, result_axpy);
    result_axpy = testvector_2;
    result_axpy.axpy(ScalarType(-0.25), countingup);
    correct_result = countingup;
    correct_result.scal(ScalarType(-0.25));
    correct_result += testvector_2;
    EXPECT_EQ(correct_result, result_axpy);
    result_axpy = testvector_4;
    result_axpy.axpy(ScalarType(-3), testvector_1);
    correct_result = testvector_1;
    correct_result.scal(ScalarType(-3));
    correct_result += testvector_4;
    EXPECT_EQ(correct_result, result_axpy);
    a = ones;
    a.axpy(ScalarType(2), testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      EXPECT_TRUE(Common::FloatCmp::eq(ScalarType(1), ones.get_entry(ii))) << "check copy-on-write";
    }
  } // void produces_correct_results() const
}; // struct VectorTest


TEST_F(VectorTest_{{T_NAME}}, fulfills_interface)
{
  this->fulfills_interface();
}
TEST_F(VectorTest_{{T_NAME}}, produces_correct_results)
{
  this->produces_correct_results();
}
{% endfor %}
