// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#include "main.hxx"

#include <type_traits>
#include <memory>

#include <dune/common/float_cmp.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/la/container/common.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/container/istl.hh>
#include <dune/stuff/la/container.hh>

#include "la_container.hh"

using namespace Dune;

static const size_t dim = 4;

typedef testing::Types<Dune::Stuff::LA::CommonDenseVector<double>
#if HAVE_EIGEN
                       ,
                       Dune::Stuff::LA::EigenDenseVector<double>, Dune::Stuff::LA::EigenMappedDenseVector<double>
#endif
#if HAVE_DUNE_ISTL
                       ,
                       Dune::Stuff::LA::IstlDenseVector<double>
#endif
                       > VectorTypes;

typedef testing::
    Types<std::pair<Dune::Stuff::LA::CommonDenseMatrix<double>, Dune::Stuff::LA::CommonDenseVector<double>>
#if HAVE_EIGEN
          ,
          std::pair<Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>, Dune::Stuff::LA::EigenDenseVector<double>>,
          std::pair<Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>,
                    Dune::Stuff::LA::EigenMappedDenseVector<double>>,
          std::pair<Dune::Stuff::LA::EigenDenseMatrix<double>, Dune::Stuff::LA::EigenDenseVector<double>>,
          std::pair<Dune::Stuff::LA::EigenDenseMatrix<double>, Dune::Stuff::LA::EigenMappedDenseVector<double>>
#endif
#if HAVE_DUNE_ISTL
          ,
          std::pair<Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>, Dune::Stuff::LA::IstlDenseVector<double>>
#endif
          > MatrixVectorCombinations;

typedef testing::Types<Dune::Stuff::LA::CommonDenseVector<double>, Dune::Stuff::LA::CommonDenseMatrix<double>
#if HAVE_EIGEN
                       ,
                       Dune::Stuff::LA::EigenDenseVector<double>, Dune::Stuff::LA::EigenRowMajorSparseMatrix<double>
#endif
#if HAVE_DUNE_ISTL
                       ,
                       Dune::Stuff::LA::IstlDenseVector<double>, Dune::Stuff::LA::IstlRowMajorSparseMatrix<double>
#endif
                       > ContainerTypes;


template <class ContainerImp>
struct ContainerTest : public ::testing::Test
{
  void fulfills_interface() const
  {
    // static tests
    typedef typename ContainerImp::Traits Traits;
    // * of the traits
    typedef typename Traits::derived_type T_derived_type;
    static_assert(std::is_same<ContainerImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    typedef typename Traits::ScalarType T_ScalarType;
    // * of the container as itself (aka the derived type)
    typedef typename ContainerImp::ScalarType D_ScalarType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // * of the container as the interface
    typedef typename Stuff::LA::ContainerInterface<Traits> InterfaceType;
    typedef typename InterfaceType::derived_type I_derived_type;
    typedef typename InterfaceType::ScalarType I_ScalarType;
    static_assert(std::is_same<ContainerImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
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
    if (!d_by_size.has_equal_shape(d_deep_copy))
      DUNE_THROW(Dune::Exception, "");
    // * of the container as the interface
    InterfaceType& i_by_size = static_cast<InterfaceType&>(d_by_size);
    ContainerImp i_deep_copy = i_by_size.copy();
    i_by_size.scal(I_ScalarType(1));
    i_by_size.axpy(I_ScalarType(1), i_deep_copy);
  } // void fulfills_interface() const
}; // struct ContainerTest


TYPED_TEST_CASE(ContainerTest, ContainerTypes);
TYPED_TEST(ContainerTest, fulfills_interface)
{
  this->fulfills_interface();
}


template <class VectorImp>
struct VectorTest : public ::testing::Test
{
  void fulfills_interface() const
  {
    // static tests
    typedef typename VectorImp::Traits Traits;
    // * of the traits
    typedef typename Traits::derived_type T_derived_type;
    static_assert(std::is_same<VectorImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    typedef typename Traits::ScalarType T_ScalarType;
    // * of the vector as itself (aka the derived type)
    typedef typename VectorImp::ScalarType D_ScalarType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // * of the vector as the interface
    typedef typename Stuff::LA::VectorInterface<Traits> InterfaceType;
    typedef typename InterfaceType::derived_type I_derived_type;
    typedef typename InterfaceType::ScalarType I_ScalarType;
    static_assert(std::is_same<VectorImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // dynamic tests
    // * of the vector as itself (aka the derived type)
    VectorImp d_by_size(dim);
    VectorImp d_by_size_and_value(dim, D_ScalarType(0));
    size_t d_size = d_by_size.size();
    if (d_size != dim)
      DUNE_THROW(Dune::Exception, d_size << " vs. " << dim);
    for (size_t ii = 0; ii < d_size; ++ii) {
      d_by_size_and_value.set_entry(ii, D_ScalarType(0.5) + D_ScalarType(ii));
      d_by_size_and_value.add_to_entry(ii, D_ScalarType(0.5) + D_ScalarType(ii));
      if (FloatCmp::ne(d_by_size_and_value.get_entry(ii), D_ScalarType(2) * D_ScalarType(ii) + D_ScalarType(1)))
        DUNE_THROW(Dune::Exception, d_by_size_and_value.get_entry(ii));
      if (FloatCmp::ne(d_by_size_and_value.get_entry(ii), d_by_size_and_value[ii]))
        DUNE_THROW(Dune::Exception, d_by_size_and_value[ii]);
    }
    size_t d_dim = d_by_size.dim();
    if (d_dim != dim)
      DUNE_THROW(Dune::Exception, d_dim << " vs. " << dim);
    bool d_almost_equal = d_by_size.almost_equal(d_by_size);
    if (!d_almost_equal)
      DUNE_THROW(Dune::Exception, "");
    d_by_size_and_value.scal(D_ScalarType(0));
    if (!d_by_size_and_value.almost_equal(d_by_size))
      DUNE_THROW(Dune::Exception, "");
    D_ScalarType d_dot = d_by_size.dot(d_by_size_and_value);
    if (!Dune::FloatCmp::eq(d_dot, D_ScalarType(0)))
      DUNE_THROW(Dune::Exception, d_dot);
    D_ScalarType d_l1_norm = d_by_size.l1_norm();
    if (!Dune::FloatCmp::eq(d_l1_norm, D_ScalarType(0)))
      DUNE_THROW(Dune::Exception, d_l1_norm);
    D_ScalarType d_l2_norm = d_by_size.l2_norm();
    if (!Dune::FloatCmp::eq(d_l2_norm, D_ScalarType(0)))
      DUNE_THROW(Dune::Exception, d_l2_norm);
    D_ScalarType d_sup_norm = d_by_size.sup_norm();
    if (!Dune::FloatCmp::eq(d_sup_norm, D_ScalarType(0)))
      DUNE_THROW(Dune::Exception, d_sup_norm);
    VectorImp d_ones(dim, D_ScalarType(1));
    std::pair<size_t, D_ScalarType> d_amax = d_ones.amax();
    if (d_amax.first != 0 || !Dune::FloatCmp::eq(d_amax.second, D_ScalarType(1)))
      DUNE_THROW(Dune::Exception, d_amax.first << ",  " << d_amax.second);
    d_ones.add(d_by_size, d_by_size_and_value);
    if (!d_by_size_and_value.almost_equal(d_ones))
      DUNE_THROW(Dune::Exception, "");
    VectorImp d_added = d_ones.add(d_by_size);
    if (!d_added.almost_equal(d_by_size_and_value))
      DUNE_THROW(Dune::Exception, "");
    d_added.iadd(d_by_size);
    if (!d_added.almost_equal(d_ones))
      DUNE_THROW(Dune::Exception, "");
    d_ones.sub(d_by_size, d_by_size_and_value);
    if (!d_by_size_and_value.almost_equal(d_ones))
      DUNE_THROW(Dune::Exception, "");
    VectorImp d_subtracted = d_ones.sub(d_by_size);
    if (!d_subtracted.almost_equal(d_by_size_and_value))
      DUNE_THROW(Dune::Exception, "");
    d_subtracted.isub(d_by_size);
    if (!d_subtracted.almost_equal(d_ones))
      DUNE_THROW(Dune::Exception, "");
    // * of the vector as the interface
    VectorImp d_by_size_2(dim);
    VectorImp d_by_size_and_value_2(dim, D_ScalarType(1));
    InterfaceType& i_by_size           = static_cast<InterfaceType&>(d_by_size_2);
    InterfaceType& i_by_size_and_value = static_cast<InterfaceType&>(d_by_size_and_value_2);
    DUNE_STUFF_SSIZE_T DUNE_UNUSED(i_dim) = i_by_size.dim();
    bool i_almost_equal = i_by_size.almost_equal(d_by_size_2);
    if (!i_almost_equal)
      DUNE_THROW(Dune::Exception, "");
    i_by_size_and_value.scal(I_ScalarType(0));
    if (!i_by_size_and_value.almost_equal(d_by_size_2))
      DUNE_THROW(Dune::Exception, "");
    I_ScalarType i_dot = i_by_size.dot(d_by_size_and_value_2);
    if (!Dune::FloatCmp::eq(i_dot, I_ScalarType(0)))
      DUNE_THROW(Dune::Exception, i_dot);
    I_ScalarType i_l1_norm = i_by_size.l1_norm();
    if (!Dune::FloatCmp::eq(i_l1_norm, I_ScalarType(0)))
      DUNE_THROW(Dune::Exception, i_l1_norm);
    I_ScalarType i_l2_norm = i_by_size.l2_norm();
    if (!Dune::FloatCmp::eq(i_l2_norm, I_ScalarType(0)))
      DUNE_THROW(Dune::Exception, i_l2_norm);
    I_ScalarType i_sup_norm = i_by_size.sup_norm();
    if (!Dune::FloatCmp::eq(i_sup_norm, I_ScalarType(0)))
      DUNE_THROW(Dune::Exception, i_sup_norm);
    VectorImp i_ones(dim, I_ScalarType(1));
    std::pair<size_t, I_ScalarType> i_amax = i_ones.amax();
    if (i_amax.first != 0 || !Dune::FloatCmp::eq(i_amax.second, I_ScalarType(1)))
      DUNE_THROW(Dune::Exception, "");
    i_ones.add(d_by_size_2, d_by_size_and_value_2);
    if (!i_by_size_and_value.almost_equal(i_ones))
      DUNE_THROW(Dune::Exception, "");
    VectorImp i_added = i_ones.add(d_by_size_2);
    if (!i_added.almost_equal(d_by_size_and_value_2))
      DUNE_THROW(Dune::Exception, "");
    i_added.iadd(d_by_size_2);
    if (!i_added.almost_equal(i_ones))
      DUNE_THROW(Dune::Exception, "");
    i_ones.sub(d_by_size_2, d_by_size_and_value_2);
    if (!i_by_size_and_value.almost_equal(i_ones))
      DUNE_THROW(Dune::Exception, "");
    VectorImp i_subtracted = i_ones.sub(d_by_size_2);
    if (!i_subtracted.almost_equal(d_by_size_and_value_2))
      DUNE_THROW(Dune::Exception, "");
    i_subtracted.isub(d_by_size_2);
    if (!i_subtracted.almost_equal(i_ones))
      DUNE_THROW(Dune::Exception, "");
  } // void fulfills_interface() const

  void produces_correct_results() const
  {
    typedef typename VectorImp::ScalarType ScalarType;

    // create test vectors
    VectorImp zeros(dim); // [0, 0, 0, 0]
    VectorImp ones(dim, ScalarType(1)); // [1, 1, 1, 1]
    VectorImp countingup(dim); //[0, 1, 2, 3]
    for (size_t ii = 0; ii < dim; ++ii)
      countingup.set_entry(ii, ScalarType(ii));
    VectorImp testvector_1(dim); // [0, -2, 2, 1]
    testvector_1[0] = ScalarType(0);
    testvector_1[1] = ScalarType(-2);
    testvector_1[2] = ScalarType(2);
    testvector_1[3] = ScalarType(1);
    VectorImp testvector_2(dim); // [0, 2, -2, 1]
    testvector_2[0] = ScalarType(0);
    testvector_2[1] = ScalarType(2);
    testvector_2[2] = ScalarType(-2);
    testvector_2[3] = ScalarType(1);
    VectorImp testvector_3(dim); // [-1, 1, -1, 1]
    testvector_3[0] = ScalarType(-1);
    testvector_3[1] = ScalarType(1);
    testvector_3[2] = ScalarType(-1);
    testvector_3[3] = ScalarType(1);
    VectorImp testvector_4(dim); // [0, 3, -2, 0]
    testvector_4[0] = ScalarType(0);
    testvector_4[1] = ScalarType(3);
    testvector_4[2] = ScalarType(-2);
    testvector_4[3] = ScalarType(0);
    VectorImp testvector_5(dim); // [1.25, 0, 2.5, -3.5]
    testvector_5[0] = ScalarType(1.25);
    testvector_5[1] = ScalarType(0);
    testvector_5[2] = ScalarType(2.5);
    testvector_5[3] = ScalarType(-3.5);

    // test amax()
    std::pair<size_t, ScalarType> amax = zeros.amax();
    if (amax.first != 0 || !Dune::FloatCmp::eq(amax.second, ScalarType(0)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = ones.amax();
    if (amax.first != 0 || !Dune::FloatCmp::eq(amax.second, ScalarType(1)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = countingup.amax();
    if (amax.first != 3 || !Dune::FloatCmp::eq(amax.second, ScalarType(3)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = testvector_1.amax();
    if (amax.first != 1 || !Dune::FloatCmp::eq(amax.second, ScalarType(2)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = testvector_2.amax();
    if (amax.first != 1 || !Dune::FloatCmp::eq(amax.second, ScalarType(2)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = testvector_3.amax();
    if (amax.first != 0 || !Dune::FloatCmp::eq(amax.second, ScalarType(1)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = testvector_4.amax();
    if (amax.first != 1 || !Dune::FloatCmp::eq(amax.second, ScalarType(3)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);
    amax = testvector_5.amax();
    if (amax.first != 3 || !Dune::FloatCmp::eq(amax.second, ScalarType(3.5)))
      DUNE_THROW(Dune::Exception, amax.first << " , " << amax.second);

    // test l1_norm()
    ScalarType l1_norm = zeros.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(0)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = ones.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(4)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = countingup.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(6)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = testvector_1.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(5)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = testvector_2.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(5)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = testvector_3.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(4)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = testvector_4.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(5)))
      DUNE_THROW(Dune::Exception, l1_norm);
    l1_norm = testvector_5.l1_norm();
    if (!Dune::FloatCmp::eq(l1_norm, ScalarType(7.25)))
      DUNE_THROW(Dune::Exception, l1_norm);

    // test l2_norm()
    ScalarType l2_norm = zeros.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(0)))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = ones.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(2)))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = countingup.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(std::sqrt(14))))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = testvector_1.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(3)))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = testvector_2.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(3)))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = testvector_3.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(2)))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = testvector_4.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(std::sqrt(13))))
      DUNE_THROW(Dune::Exception, l2_norm);
    l2_norm = testvector_5.l2_norm();
    if (!Dune::FloatCmp::eq(l2_norm, ScalarType(std::sqrt(20.0625))))
      DUNE_THROW(Dune::Exception, l2_norm);

    // test sup_norm()
    ScalarType sup_norm = zeros.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(0)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = ones.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(1)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = countingup.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(3)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = testvector_1.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(2)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = testvector_2.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(2)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = testvector_3.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(1)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = testvector_4.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(3)))
      DUNE_THROW(Dune::Exception, sup_norm);
    sup_norm = testvector_5.sup_norm();
    if (!Dune::FloatCmp::eq(sup_norm, ScalarType(3.5)))
      DUNE_THROW(Dune::Exception, sup_norm);

    // test dot(), operator*
    ScalarType dot            = ones.dot(zeros);
    ScalarType dot_operator   = ones * zeros;
    ScalarType dot2           = zeros.dot(ones);
    ScalarType dot_operator_2 = zeros * ones;
    if (!Dune::FloatCmp::eq(dot, ScalarType(0)) || !Dune::FloatCmp::eq(dot, dot2)
        || !Dune::FloatCmp::eq(dot_operator, dot_operator_2)
        || !Dune::FloatCmp::eq(dot, dot_operator))
      DUNE_THROW(Dune::Exception, dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2);
    dot          = ones.dot(ones);
    dot_operator = ones * ones;
    if (!Dune::FloatCmp::eq(dot, ScalarType(4)) || !Dune::FloatCmp::eq(dot_operator, ScalarType(4)))
      DUNE_THROW(Dune::Exception, dot << ", " << dot_operator);
    dot            = ones.dot(testvector_3);
    dot_operator   = ones * testvector_3;
    dot2           = testvector_3.dot(ones);
    dot_operator_2 = testvector_3 * ones;
    if (!Dune::FloatCmp::eq(dot, ScalarType(0)) || !Dune::FloatCmp::eq(dot, dot2)
        || !Dune::FloatCmp::eq(dot_operator, dot_operator_2)
        || !Dune::FloatCmp::eq(dot, dot_operator))
      DUNE_THROW(Dune::Exception, dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2);
    dot            = countingup.dot(testvector_5);
    dot_operator   = countingup * testvector_5;
    dot2           = testvector_5.dot(countingup);
    dot_operator_2 = testvector_5 * countingup;
    if (!Dune::FloatCmp::eq(dot, ScalarType(-5.5)) || !Dune::FloatCmp::eq(dot, dot2)
        || !Dune::FloatCmp::eq(dot_operator, dot_operator_2)
        || !Dune::FloatCmp::eq(dot, dot_operator))
      DUNE_THROW(Dune::Exception, dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2);
    dot            = testvector_3.dot(testvector_5);
    dot_operator   = testvector_3 * testvector_5;
    dot2           = testvector_5.dot(testvector_3);
    dot_operator_2 = testvector_5 * testvector_3;
    if (!Dune::FloatCmp::eq(dot, ScalarType(-7.25)) || !Dune::FloatCmp::eq(dot, dot2)
        || !Dune::FloatCmp::eq(dot_operator, dot_operator_2)
        || !Dune::FloatCmp::eq(dot, dot_operator))
      DUNE_THROW(Dune::Exception, dot << ", " << dot2 << ", " << dot_operator << ", " << dot_operator_2);

    // test operator==
    if (!(zeros == zeros))
      DUNE_THROW(Dune::Exception, "");
    if (!(ones == ones))
      DUNE_THROW(Dune::Exception, "");
    if (!(countingup == countingup))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_1 == testvector_1))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_2 == testvector_2))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_3 == testvector_3))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_4 == testvector_4))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_5 == testvector_5))
      DUNE_THROW(Dune::Exception, "");

    // test operator!=
    if (!(zeros != ones))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_1 != testvector_2))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_3 != testvector_4))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_3 != ones))
      DUNE_THROW(Dune::Exception, "");
    if (!(testvector_4 != testvector_5))
      DUNE_THROW(Dune::Exception, "");

    // test scal, operator*=
    VectorImp scaled             = zeros;
    VectorImp scaled_by_operator = zeros;
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    if (scaled != zeros || scaled_by_operator != zeros)
      DUNE_THROW(Dune::Exception, "");
    scaled             = ones;
    scaled_by_operator = ones;
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    if (scaled != ones || scaled_by_operator != ones)
      DUNE_THROW(Dune::Exception, "");
    scaled.scal(ScalarType(0));
    scaled_by_operator *= ScalarType(0);
    if (scaled != zeros || scaled_by_operator != zeros)
      DUNE_THROW(Dune::Exception, "");
    scaled = testvector_1;
    scaled.scal(ScalarType(2));
    scaled_by_operator = testvector_1;
    scaled_by_operator *= ScalarType(2);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(scaled[ii], ScalarType(2) * testvector_1[ii])
          || !Dune::FloatCmp::eq(scaled_by_operator[ii], ScalarType(2) * testvector_1[ii]))
        DUNE_THROW(Dune::Exception, scaled[ii] << " vs. " << scaled_by_operator[ii] << " vs." << testvector_1[ii]);
    }
    scaled = testvector_3;
    scaled.scal(ScalarType(-2));
    scaled_by_operator = testvector_3;
    scaled_by_operator *= ScalarType(-2);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(scaled[ii], ScalarType(-2) * testvector_3[ii])
          || !Dune::FloatCmp::eq(scaled_by_operator[ii], ScalarType(-2) * testvector_3[ii]))
        DUNE_THROW(Dune::Exception, scaled[ii] << " vs. " << scaled_by_operator[ii] << " vs." << testvector_3[ii]);
    }
    scaled = countingup;
    scaled.scal(ScalarType(2.2));
    scaled_by_operator = countingup;
    scaled_by_operator *= ScalarType(2.2);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(scaled[ii], ScalarType(2.2) * countingup[ii])
          || !Dune::FloatCmp::eq(scaled_by_operator[ii], ScalarType(2.2) * countingup[ii]))
        DUNE_THROW(Dune::Exception, scaled[ii] << " vs. " << scaled_by_operator[ii] << " vs." << countingup[ii]);
    }
    scaled = testvector_5;
    scaled.scal(ScalarType(-3.75));
    scaled_by_operator = testvector_5;
    scaled_by_operator *= ScalarType(-3.75);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(scaled[ii], ScalarType(-3.75) * testvector_5[ii])
          || !Dune::FloatCmp::eq(scaled_by_operator[ii], ScalarType(-3.75) * testvector_5[ii]))
        DUNE_THROW(Dune::Exception, scaled[ii] << " vs. " << scaled_by_operator[ii] << " vs." << testvector_5[ii]);
    }
    VectorImp a = ones;
    a.scal(ScalarType(0));
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }

    // test operator+, operator+=, add, iadd
    VectorImp sum_operator_plus  = zeros + ones;
    VectorImp sum_operator_iplus = zeros;
    sum_operator_iplus += ones;
    VectorImp sum_add_1 = zeros.add(ones);
    VectorImp sum_add_2(dim);
    zeros.add(ones, sum_add_2);
    VectorImp sum_iadd = zeros;
    sum_iadd.iadd(ones);
    VectorImp sum_correct(dim);
    sum_correct[0] = ScalarType(1);
    sum_correct[1] = ScalarType(1);
    sum_correct[2] = ScalarType(1);
    sum_correct[3] = ScalarType(1);
    if (sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
        || sum_add_1 != sum_operator_plus
        || sum_add_1 != sum_correct)
      DUNE_THROW(Dune::Exception,
                 sum_add_1 << ", " << sum_add_2 << ", " << sum_operator_iplus << ", " << sum_operator_plus << ", "
                           << sum_iadd
                           << ", "
                           << sum_correct);

    sum_operator_plus  = countingup + testvector_1;
    sum_operator_iplus = countingup;
    sum_operator_iplus += testvector_1;
    sum_add_1 = countingup.add(testvector_1);
    countingup.add(testvector_1, sum_add_2);
    sum_iadd = countingup;
    sum_iadd.iadd(testvector_1);
    sum_correct[0] = ScalarType(0);
    sum_correct[1] = ScalarType(-1);
    sum_correct[2] = ScalarType(4);
    sum_correct[3] = ScalarType(4);
    if (sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
        || sum_add_1 != sum_operator_plus
        || sum_add_1 != sum_correct)
      DUNE_THROW(Dune::Exception,
                 sum_add_1 << ", " << sum_add_2 << ", " << sum_operator_iplus << ", " << sum_operator_plus << ", "
                           << sum_iadd
                           << ", "
                           << sum_correct);

    sum_operator_plus  = testvector_3 + testvector_5;
    sum_operator_iplus = testvector_3;
    sum_operator_iplus += testvector_5;
    sum_add_1 = testvector_3.add(testvector_5);
    testvector_3.add(testvector_5, sum_add_2);
    sum_iadd = testvector_3;
    sum_iadd.iadd(testvector_5);
    sum_correct[0] = ScalarType(0.25);
    sum_correct[1] = ScalarType(1);
    sum_correct[2] = ScalarType(1.5);
    sum_correct[3] = ScalarType(-2.5);
    if (sum_operator_iplus != sum_operator_plus || sum_add_1 != sum_add_2 || sum_iadd != sum_add_1
        || sum_add_1 != sum_operator_plus
        || sum_add_1 != sum_correct)
      DUNE_THROW(Dune::Exception,
                 sum_add_1 << ", " << sum_add_2 << ", " << sum_operator_iplus << ", " << sum_operator_plus << ", "
                           << sum_iadd
                           << ", "
                           << sum_correct);

    a = ones;
    a += testvector_3;
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }
    a = ones;
    a.iadd(testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }
    a = ones;
    ones.add(testvector_3, a);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }

    // test operator-, operator-=, sub, isub
    VectorImp diff_operator_minus  = zeros - ones;
    VectorImp diff_operator_iminus = zeros;
    diff_operator_iminus -= ones;
    VectorImp diff_sub_1 = zeros.sub(ones);
    VectorImp diff_sub_2(dim);
    zeros.sub(ones, diff_sub_2);
    VectorImp diff_isub = zeros;
    diff_isub.isub(ones);
    VectorImp diff_correct(dim);
    diff_correct[0] = ScalarType(-1);
    diff_correct[1] = ScalarType(-1);
    diff_correct[2] = ScalarType(-1);
    diff_correct[3] = ScalarType(-1);
    if (diff_operator_iminus != diff_operator_minus || diff_sub_1 != diff_sub_2 || diff_isub != diff_sub_1
        || diff_sub_1 != diff_operator_minus
        || diff_sub_1 != diff_correct)
      DUNE_THROW(Dune::Exception,
                 diff_sub_1 << ", " << diff_sub_2 << ", " << diff_operator_iminus << ", " << diff_operator_minus << ", "
                            << diff_isub
                            << ", "
                            << diff_correct);

    diff_operator_minus  = testvector_1 - testvector_4;
    diff_operator_iminus = testvector_1;
    diff_operator_iminus -= testvector_4;
    diff_sub_1 = testvector_1.sub(testvector_4);
    testvector_1.sub(testvector_4, diff_sub_2);
    diff_isub = testvector_1;
    diff_isub.isub(testvector_4);
    diff_correct[0] = ScalarType(-0);
    diff_correct[1] = ScalarType(-5);
    diff_correct[2] = ScalarType(4);
    diff_correct[3] = ScalarType(1);
    if (diff_operator_iminus != diff_operator_minus || diff_sub_1 != diff_sub_2 || diff_isub != diff_sub_1
        || diff_sub_1 != diff_operator_minus
        || diff_sub_1 != diff_correct)
      DUNE_THROW(Dune::Exception,
                 diff_sub_1 << ", " << diff_sub_2 << ", " << diff_operator_iminus << ", " << diff_operator_minus << ", "
                            << diff_isub
                            << ", "
                            << diff_correct);

    diff_operator_minus  = testvector_5 - testvector_2;
    diff_operator_iminus = testvector_5;
    diff_operator_iminus -= testvector_2;
    diff_sub_1 = testvector_5.sub(testvector_2);
    testvector_5.sub(testvector_2, diff_sub_2);
    diff_isub = testvector_5;
    diff_isub.isub(testvector_2);
    diff_correct[0] = ScalarType(1.25);
    diff_correct[1] = ScalarType(-2);
    diff_correct[2] = ScalarType(4.5);
    diff_correct[3] = ScalarType(-4.5);
    if (diff_operator_iminus != diff_operator_minus || diff_sub_1 != diff_sub_2 || diff_isub != diff_sub_1
        || diff_sub_1 != diff_operator_minus
        || diff_sub_1 != diff_correct)
      DUNE_THROW(Dune::Exception,
                 diff_sub_1 << ", " << diff_sub_2 << ", " << diff_operator_iminus << ", " << diff_operator_minus << ", "
                            << diff_isub
                            << ", "
                            << diff_correct);

    a = ones;
    a -= testvector_3;
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }
    a = ones;
    a.isub(testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }
    a = ones;
    ones.sub(testvector_3, a);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }

    // test axpy
    VectorImp result_axpy = zeros;
    result_axpy.axpy(ScalarType(1), ones);
    if (result_axpy != ones)
      DUNE_THROW(Dune::Exception, result_axpy << " vs. " << ones);
    result_axpy = testvector_5;
    result_axpy.axpy(ScalarType(2.75), testvector_3);
    VectorImp correct_result = testvector_3;
    correct_result.scal(ScalarType(2.75));
    correct_result += testvector_5;
    if (result_axpy != correct_result)
      DUNE_THROW(Dune::Exception, result_axpy << " vs. " << correct_result);
    result_axpy = testvector_2;
    result_axpy.axpy(ScalarType(-0.25), countingup);
    correct_result = countingup;
    correct_result.scal(ScalarType(-0.25));
    correct_result += testvector_2;
    if (result_axpy != correct_result)
      DUNE_THROW(Dune::Exception, result_axpy << " vs. " << correct_result);
    result_axpy = testvector_4;
    result_axpy.axpy(ScalarType(-3), testvector_1);
    correct_result = testvector_1;
    correct_result.scal(ScalarType(-3));
    correct_result += testvector_4;
    if (result_axpy != correct_result)
      DUNE_THROW(Dune::Exception, result_axpy << " vs. " << correct_result);
    a = ones;
    a.axpy(ScalarType(2), testvector_3);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }
  } // void produces_correct_results() const
}; // struct VectorTest


TYPED_TEST_CASE(VectorTest, VectorTypes);
TYPED_TEST(VectorTest, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST_CASE(VectorTest, VectorTypes);
TYPED_TEST(VectorTest, produces_correct_results)
{
  this->produces_correct_results();
}

template <class MatrixVectorCombination>
struct MatrixTest : public ::testing::Test
{
  typedef typename MatrixVectorCombination::first_type MatrixImp;
  typedef typename MatrixVectorCombination::second_type VectorImp;
  typedef typename Dune::Stuff::LA::SparsityPatternDefault PatternType;

  void fulfills_interface() const
  {
    // static tests
    typedef typename MatrixImp::Traits Traits;
    // * of the traits
    typedef typename Traits::derived_type T_derived_type;
    static_assert(std::is_same<MatrixImp, T_derived_type>::value, "derived_type has to be the correct Type!");
    typedef typename Traits::ScalarType T_ScalarType;
    // * of the matrix as itself (aka the derived type)
    typedef typename MatrixImp::ScalarType D_ScalarType;
    static_assert(std::is_same<T_ScalarType, D_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // * of the matrix as the interface
    typedef typename Stuff::LA::MatrixInterface<Traits> InterfaceType;
    typedef typename InterfaceType::derived_type I_derived_type;
    typedef typename InterfaceType::ScalarType I_ScalarType;
    static_assert(std::is_same<MatrixImp, I_derived_type>::value, "derived_type has to be the correct Type!");
    static_assert(std::is_same<T_ScalarType, I_ScalarType>::value,
                  "ScalarType of derived_type has to be the correct Type!");
    // dynamic tests
    // * of the matrix as itself (aka the derived type)
    MatrixImp d_by_size(dim, dim);
    PatternType pattern(dim);
    for (size_t ii = 0; ii < dim; ++ii) {
      for (size_t jj = 0; jj < dim; ++jj)
        pattern.inner(ii).push_back(jj);
    }
    MatrixImp d_by_size_and_pattern(dim, dim, pattern);
    size_t d_rows = d_by_size.rows();
    if (d_rows != dim)
      DUNE_THROW(Dune::Exception, d_rows << " vs. " << dim);
    size_t d_cols = d_by_size.cols();
    if (d_cols != dim)
      DUNE_THROW(Dune::Exception, d_cols << " vs. " << dim);
    d_rows = d_by_size_and_pattern.rows();
    if (d_rows != dim)
      DUNE_THROW(Dune::Exception, d_rows << " vs. " << dim);
    d_cols = d_by_size_and_pattern.cols();
    if (d_cols != dim)
      DUNE_THROW(Dune::Exception, d_cols << " vs. " << dim);
    VectorImp zeros(dim);
    VectorImp ones(dim, D_ScalarType(1));
    VectorImp result(dim);
    d_by_size_and_pattern.mv(ones, result);
    if (result != zeros)
      DUNE_THROW(Dune::Exception, "check mv");
    for (size_t ii = 0; ii < d_rows; ++ii) {
      d_by_size_and_pattern.unit_row(ii);
      if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, ii), D_ScalarType(1)))
        DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, ii) << " vs " << D_ScalarType(1));
      for (size_t jj = 0; jj < ii; ++jj) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
      for (size_t jj = ii + 1; jj < d_cols; ++jj) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
    }
    for (size_t ii = 0; ii < d_rows; ++ii) {
      d_by_size_and_pattern.clear_row(ii);
      for (size_t jj = 0; jj < d_cols; ++jj) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
    }
    for (size_t jj = 0; jj < d_cols; ++jj) {
      d_by_size_and_pattern.unit_col(jj);
      if (FloatCmp::ne(d_by_size_and_pattern.get_entry(jj, jj), D_ScalarType(1)))
        DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(jj, jj) << " vs " << D_ScalarType(1));
      for (size_t ii = 0; ii < jj; ++ii) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
      for (size_t ii = jj + 1; ii < d_rows; ++ii) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
    }
    for (size_t jj = 0; jj < d_cols; ++jj) {
      d_by_size_and_pattern.clear_col(jj);
      for (size_t ii = 0; ii < d_rows; ++ii) {
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj), D_ScalarType(0)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj) << " vs " << D_ScalarType(0));
      }
    }
    for (size_t ii = 0; ii < d_rows; ++ii) {
      for (size_t jj = 0; jj < d_cols; ++jj) {
        d_by_size_and_pattern.set_entry(ii, jj, D_ScalarType(0.5) + D_ScalarType(ii) + D_ScalarType(jj));
        d_by_size_and_pattern.add_to_entry(ii, jj, D_ScalarType(0.5) + D_ScalarType(ii) + D_ScalarType(jj));
        if (FloatCmp::ne(d_by_size_and_pattern.get_entry(ii, jj),
                         D_ScalarType(2) * D_ScalarType(ii) + D_ScalarType(2) * D_ScalarType(jj) + D_ScalarType(1)))
          DUNE_THROW(Dune::Exception, d_by_size_and_pattern.get_entry(ii, jj));
      }
    }
  } // void fulfills_interface() const

  void produces_correct_results() const
  {
    typedef typename MatrixImp::ScalarType ScalarType;

    // create test patterns
    PatternType dense_pattern(dim);
    for (size_t ii = 0; ii < dim; ++ii) {
      for (size_t jj = 0; jj < dim; ++jj)
        dense_pattern.inner(ii).push_back(jj);
    }
    PatternType sparse_pattern(dim);
    sparse_pattern.inner(0).push_back(2); //|-, -, x, -|
    sparse_pattern.inner(1).push_back(0); //|x, x, -, -|
    sparse_pattern.inner(1).push_back(1); //|-, -, -, x|
    sparse_pattern.inner(2).push_back(3); //|-, -, -, -|

    // create test matrizes
    MatrixImp matrix_zeros_dense(dim, dim, dense_pattern); // |0, 0, 0, 0|
    for (size_t ii = 0; ii < dim; ++ii) { // |0, 0, 0, 0|
      for (size_t jj = 0; jj < dim; ++jj) // |0, 0, 0, 0|
        matrix_zeros_dense.set_entry(ii, jj, ScalarType(0)); // |0, 0, 0, 0|
    }
    MatrixImp matrix_zeros_sparse(dim, dim, sparse_pattern); //|-, -, 0, -|
    matrix_zeros_sparse.set_entry(0, 2, ScalarType(0)); //|0, 0, -, -|
    matrix_zeros_sparse.set_entry(1, 0, ScalarType(0)); //|-, -, -, 0|
    matrix_zeros_sparse.set_entry(1, 1, ScalarType(0)); //|-, -, -, -|
    matrix_zeros_sparse.set_entry(2, 3, ScalarType(0));
    MatrixImp matrix_ones(dim, dim, dense_pattern); // |1, 1, 1, 1|
    for (size_t ii = 0; ii < dim; ++ii) { // |1, 1, 1, 1|
      for (size_t jj = 0; jj < dim; ++jj) // |1, 1, 1, 1|
        matrix_ones.set_entry(ii, jj, ScalarType(1)); // |1, 1, 1, 1|
    }
    MatrixImp testmatrix_1(dim, dim, dense_pattern); // |0, 1, 2, 3|
    for (size_t ii = 0; ii < dim; ++ii) { // |1, 2, 3, 4|
      for (size_t jj = 0; jj < dim; ++jj) // |2, 3, 4, 5|
        testmatrix_1.set_entry(ii, jj, ScalarType(ii) + ScalarType(jj)); // |3, 4, 5, 6|
    }
    MatrixImp testmatrix_2(dim, dim, dense_pattern); // | 0.0, 1.5, 3.0, 4.5|
    for (size_t ii = 0; ii < dim; ++ii) { // |-0.5, 1.0, 2.5, 4.0|
      for (size_t jj = 0; jj < dim; ++jj) // |-1.0, 0.5, 2.0, 3.5|
        testmatrix_2.set_entry(
            ii, jj, ScalarType(-0.5) * ScalarType(ii) + ScalarType(1.5) * ScalarType(jj)); // |-1.5, 0.0, 1.5, 3.0|
    }
    MatrixImp testmatrix_sparse(dim, dim, sparse_pattern); //|-,   -, 0.5,    -|
    testmatrix_sparse.set_entry(0, 2, ScalarType(0.5)); //|1, 1.5,   -,    -|
    testmatrix_sparse.set_entry(1, 0, ScalarType(1)); //|-,   -,   -, -0.5|
    testmatrix_sparse.set_entry(1, 1, ScalarType(1.5)); //|-,   -,   -,    -|
    testmatrix_sparse.set_entry(2, 3, ScalarType(-0.5));

    // create test vectors
    VectorImp vector_zeros(dim); // [0, 0, 0, 0]
    VectorImp vector_ones(dim, ScalarType(1)); // [1, 1, 1, 1]
    VectorImp vector_countingup(dim); //[0, 1, 2, 3]
    for (size_t ii = 0; ii < dim; ++ii)
      vector_countingup.set_entry(ii, ScalarType(ii));
    VectorImp testvector_1(dim); // [0, -2, 2, 1]
    testvector_1[0] = ScalarType(0);
    testvector_1[1] = ScalarType(-2);
    testvector_1[2] = ScalarType(2);
    testvector_1[3] = ScalarType(1);
    VectorImp testvector_3(dim); // [-1, 1, -1, 1]
    testvector_3[0] = ScalarType(-1);
    testvector_3[1] = ScalarType(1);
    testvector_3[2] = ScalarType(-1);
    testvector_3[3] = ScalarType(1);
    VectorImp testvector_5(dim); // [1.25, 0, 2.5, -3.5]
    testvector_5[0] = ScalarType(1.25);
    testvector_5[1] = ScalarType(0);
    testvector_5[2] = ScalarType(2.5);
    testvector_5[3] = ScalarType(-3.5);

    // test mv
    VectorImp result_mv_1(dim);
    matrix_zeros_dense.mv(vector_zeros, result_mv_1);
    VectorImp result_mv_2(dim);
    matrix_zeros_sparse.mv(vector_zeros, result_mv_2);
    if (result_mv_1 != vector_zeros || result_mv_2 != vector_zeros)
      DUNE_THROW(Dune::Exception, result_mv_1 << " vs " << result_mv_2 << " vs. " << vector_zeros);
    testmatrix_sparse.mv(testvector_5, result_mv_1);
    if (result_mv_1[0] != ScalarType(1.25) || result_mv_1[1] != ScalarType(1.25) || result_mv_1[2] != ScalarType(1.75)
        || result_mv_1[3] != ScalarType(0))
      DUNE_THROW(Dune::Exception, result_mv_1);
    testmatrix_2.mv(testvector_3, result_mv_1);
    result_mv_2 = vector_ones;
    result_mv_2.scal(ScalarType(3));
    if (result_mv_1 != result_mv_2)
      DUNE_THROW(Dune::Exception, result_mv_1 << " vs. " << result_mv_2);
    testmatrix_1.mv(testvector_1, result_mv_1);
    if (result_mv_1[0] != ScalarType(5) || result_mv_1[1] != ScalarType(6) || result_mv_1[2] != ScalarType(7)
        || result_mv_1[3] != ScalarType(8))
      DUNE_THROW(Dune::Exception, result_mv_1);
    testmatrix_sparse.mv(vector_ones, result_mv_1);
    if (result_mv_1[0] != ScalarType(0.5) || result_mv_1[1] != ScalarType(2.5) || result_mv_1[2] != ScalarType(-0.5)
        || result_mv_1[3] != ScalarType(0))
      DUNE_THROW(Dune::Exception, result_mv_1);
    VectorImp a = vector_ones;
    matrix_zeros_dense.mv(vector_zeros, a);
    for (size_t ii = 0; ii < dim; ++ii) {
      if (!Dune::FloatCmp::eq(vector_ones[ii], ScalarType(1)))
        DUNE_THROW(Dune::Exception, "check copy-on-write");
    }

    // test scal, operator*
    MatrixImp scaled             = matrix_zeros_dense;
    MatrixImp scaled_by_operator = matrix_zeros_dense;
    size_t rows                  = scaled.rows();
    size_t cols = scaled.cols();
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(scaled.get_entry(ii, jj), ScalarType(0))
            || FloatCmp::ne(scaled_by_operator.get_entry(ii, jj), ScalarType(0)))
          DUNE_THROW(Dune::Exception, "");
      }
    }
    scaled             = matrix_zeros_sparse;
    scaled_by_operator = matrix_zeros_sparse;
    scaled.scal(ScalarType(1));
    scaled_by_operator *= ScalarType(1);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(scaled.get_entry(ii, jj), ScalarType(0))
            || FloatCmp::ne(scaled_by_operator.get_entry(ii, jj), ScalarType(0)))
          DUNE_THROW(Dune::Exception, "");
      }
    }
    scaled             = matrix_ones;
    scaled_by_operator = matrix_ones;
    scaled.scal(ScalarType(0.5));
    scaled_by_operator *= ScalarType(0.5);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(scaled.get_entry(ii, jj), ScalarType(0.5))
            || FloatCmp::ne(scaled_by_operator.get_entry(ii, jj), ScalarType(0.5)))
          DUNE_THROW(Dune::Exception, "");
      }
    }
    scaled             = testmatrix_sparse;
    scaled_by_operator = testmatrix_sparse;
    scaled.scal(ScalarType(-1.25));
    scaled_by_operator *= ScalarType(-1.25);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(scaled.get_entry(ii, jj), testmatrix_sparse.get_entry(ii, jj) * ScalarType(-1.25))
            || FloatCmp::ne(scaled_by_operator.get_entry(ii, jj),
                            testmatrix_sparse.get_entry(ii, jj) * ScalarType(-1.25)))
          DUNE_THROW(Dune::Exception, "");
      }
    }
    scaled             = testmatrix_1;
    scaled_by_operator = testmatrix_1;
    scaled.scal(ScalarType(10));
    scaled_by_operator *= ScalarType(10);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(scaled.get_entry(ii, jj), testmatrix_1.get_entry(ii, jj) * ScalarType(10))
            || FloatCmp::ne(scaled_by_operator.get_entry(ii, jj), testmatrix_1.get_entry(ii, jj) * ScalarType(10)))
          DUNE_THROW(Dune::Exception, "");
      }
    }
    MatrixImp b = matrix_ones;
    b.scal(ScalarType(0));
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(matrix_ones.get_entry(ii, jj), ScalarType(1)))
          DUNE_THROW(Dune::Exception, "check copy-on-write");
      }
    }

    // test axpy
    MatrixImp result_axpy = matrix_zeros_dense;
    result_axpy.axpy(ScalarType(1.5), matrix_ones);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(result_axpy.get_entry(ii, jj), ScalarType(1.5)))
          DUNE_THROW(Dune::Exception, result_axpy.get_entry(ii, jj) << " vs. " << ScalarType(1.5));
      }
    }
    result_axpy = matrix_zeros_sparse;
    result_axpy.axpy(ScalarType(-1.5), matrix_zeros_sparse);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(result_axpy.get_entry(ii, jj), ScalarType(0)))
          DUNE_THROW(Dune::Exception, result_axpy.get_entry(ii, jj) << " vs. " << ScalarType(0));
      }
    }
    result_axpy = testmatrix_sparse;
    result_axpy.axpy(ScalarType(-0.5), testmatrix_sparse);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(result_axpy.get_entry(ii, jj), ScalarType(0.5) * testmatrix_sparse.get_entry(ii, jj)))
          DUNE_THROW(Dune::Exception,
                     result_axpy.get_entry(ii, jj) << " vs. " << ScalarType(0.5) * testmatrix_sparse.get_entry(ii, jj));
      }
    }
    result_axpy = testmatrix_1;
    result_axpy.axpy(ScalarType(2), testmatrix_2);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(result_axpy.get_entry(ii, jj),
                         ScalarType(2) * testmatrix_2.get_entry(ii, jj) + testmatrix_1.get_entry(ii, jj)))
          DUNE_THROW(Dune::Exception,
                     result_axpy.get_entry(ii, jj)
                         << " vs. "
                         << ScalarType(2) * testmatrix_2.get_entry(ii, jj) + testmatrix_1.get_entry(ii, jj));
      }
    }
    b = matrix_zeros_dense;
    b.axpy(ScalarType(1), matrix_ones);
    for (size_t ii = 0; ii < rows; ++ii) {
      for (size_t jj = 0; jj < cols; ++jj) {
        if (FloatCmp::ne(matrix_zeros_dense.get_entry(ii, jj), ScalarType(0)))
          DUNE_THROW(Dune::Exception, "check copy-on-write");
      }
    }
  } // void produces_correct_results() const
}; // struct MatrixTest

TYPED_TEST_CASE(MatrixTest, MatrixVectorCombinations);
TYPED_TEST(MatrixTest, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST_CASE(MatrixTest, MatrixVectorCombinations);
TYPED_TEST(MatrixTest, produces_correct_results)
{
  this->produces_correct_results();
}
