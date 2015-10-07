// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_TEST_LA_CONTAINER_HH
#define DUNE_STUFF_TEST_LA_CONTAINER_HH

#include <dune/stuff/la/container.hh>

template <class ContainerImp>
class ContainerFactory
{
public:
  static ContainerImp create(const size_t /*size*/)
  {
    static_assert(Dune::AlwaysFalse<ContainerImp>::value, "Please specialize this class for this ContainerImp!");
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::CommonDenseVector<S>>
{
public:
  static Dune::Stuff::LA::CommonDenseVector<S> create(const size_t size)
  {
    return Dune::Stuff::LA::CommonDenseVector<S>(size, S(1));
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::CommonDenseMatrix<S>>
{
public:
  static Dune::Stuff::LA::CommonDenseMatrix<S> create(const size_t size)
  {
    Dune::Stuff::LA::CommonDenseMatrix<S> matrix(size, size);
    for (size_t ii = 0; ii < size; ++ii)
      matrix.unit_row(ii);
    return matrix;
  }
};

#if HAVE_DUNE_ISTL
template <class S>
class ContainerFactory<Dune::Stuff::LA::IstlDenseVector<S>>
{
public:
  static Dune::Stuff::LA::IstlDenseVector<S> create(const size_t size)
  {
    return Dune::Stuff::LA::IstlDenseVector<S>(size, S(1));
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::IstlRowMajorSparseMatrix<S>>
{
public:
  static Dune::Stuff::LA::IstlRowMajorSparseMatrix<S> create(const size_t size)
  {
    Dune::Stuff::LA::SparsityPatternDefault pattern(size);
    for (size_t ii = 0; ii < size; ++ii)
      pattern.inner(ii).push_back(ii);
    Dune::Stuff::LA::IstlRowMajorSparseMatrix<S> matrix(size, size, pattern);
    for (size_t ii = 0; ii < size; ++ii)
      matrix.unit_row(ii);
    return matrix;
  }
};
#endif // HAVE_DUNE_ISTL

#if HAVE_EIGEN
template <class S>
class ContainerFactory<Dune::Stuff::LA::EigenDenseVector<S>>
{
public:
  static Dune::Stuff::LA::EigenDenseVector<S> create(const size_t size)
  {
    return Dune::Stuff::LA::EigenDenseVector<S>(size, S(1));
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::EigenMappedDenseVector<S>>
{
public:
  static Dune::Stuff::LA::EigenMappedDenseVector<S> create(const size_t size)
  {
    return Dune::Stuff::LA::EigenMappedDenseVector<S>(size, S(1));
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::EigenDenseMatrix<S>>
{
public:
  static Dune::Stuff::LA::EigenDenseMatrix<S> create(const size_t size)
  {
    Dune::Stuff::LA::EigenDenseMatrix<S> matrix(size, size);
    for (size_t ii = 0; ii < size; ++ii)
      matrix.unit_row(ii);
    return matrix;
  }
};

template <class S>
class ContainerFactory<Dune::Stuff::LA::EigenRowMajorSparseMatrix<S>>
{
public:
  static Dune::Stuff::LA::EigenRowMajorSparseMatrix<S> create(const size_t size)
  {
    Dune::Stuff::LA::SparsityPatternDefault pattern(size);
    for (size_t ii = 0; ii < size; ++ii)
      pattern.inner(ii).push_back(ii);
    Dune::Stuff::LA::EigenRowMajorSparseMatrix<S> matrix(size, size, pattern);
    for (size_t ii = 0; ii < size; ++ii)
      matrix.unit_row(ii);
    return matrix;
  }
};
#endif // HAVE_EIGEN

#endif // DUNE_STUFF_TEST_LA_CONTAINER_HH
