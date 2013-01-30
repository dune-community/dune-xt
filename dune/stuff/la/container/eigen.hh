#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#if HAVE_EIGEN

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <dune/common/shared_ptr.hh>

#include "interface.hh"
#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {


template <class Traits>
class EigenInterface : public Interface<Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
}; // class EigenInterface


template <class Traits>
class EigenMatrixInterface : public EigenInterface<Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
}; // class EigenMatrixInterface


template <class Traits>
class EigenVectorInterface
{
public:
  typedef typename Traits::derived_type derived_type;
}; // class EigenVectorInterface


template <class ElementType>
class EigenRowMajorSparseMatrix;


template <class ElementImp>
class EigenDenseMatrix;


template <class ElementImp>
class EigenDenseVector;


template <class ElementImp = double>
class EigenRowMajorSparseMatrixTraits
{
public:
  typedef ElementImp ElementType;
  typedef EigenRowMajorSparseMatrix<ElementType> derived_type;
  //  typedef derived_type EigenMatrix_derived_type;
  typedef typename ::Eigen::SparseMatrix<ElementType, ::Eigen::RowMajor> BackendType;
  typedef typename BackendType::Index size_type;
}; // class RowMajorSparseMatrixTraits


template <class ElementImp = double>
class EigenRowMajorSparseMatrix : public MatrixInterface<EigenRowMajorSparseMatrixTraits<ElementImp>>,
                                  public EigenMatrixInterface<EigenRowMajorSparseMatrixTraits<ElementImp>>
{
public:
  typedef EigenRowMajorSparseMatrix<ElementImp> ThisType;
  typedef EigenRowMajorSparseMatrixTraits<ElementImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ElementType ElementType;
  typedef typename Traits::size_type size_type;

  EigenRowMajorSparseMatrix()
  {
  }

  EigenRowMajorSparseMatrix(const ThisType& _other)
    : eigenMatrix_(_other.backend())
  {
  }

  EigenRowMajorSparseMatrix(const BackendType& _otherEigenMatrix)
    : eigenMatrix_(_otherEigenMatrix)
  {
  }

  ThisType& operator=(const ThisType& _other)
  {
    eigenMatrix_ = _other.backend();
    return *this;
  }

  ThisType& operator=(const BackendType& _otherEigenMatrix)
  {
    eigenMatrix_ = _otherEigenMatrix;
    return *this;
  }

  EigenRowMajorSparseMatrix(const size_type _rows, const size_type _cols,
                            const Dune::Stuff::LA::Container::SparsityPatternDefault& _pattern)
    : eigenMatrix_(_rows, _cols)
  {
    assert(size_type(_pattern.size()) == _rows && "Given pattern too short!");
    typedef Dune::Stuff::LA::Container::SparsityPatternDefault PatternType;
    typedef PatternType::SetType ColumnsType;
    for (size_type row = 0; row < size_type(_pattern.size()); ++row) {
      eigenMatrix_.startVec(row);
      const ColumnsType& columns = _pattern.set(row);
      for (typename ColumnsType::const_iterator columnIt = columns.begin(); columnIt != columns.end(); ++columnIt) {
        const size_type column = *columnIt;
        eigenMatrix_.insertBackByOuterInner(row, column);
      }
      // create diagonal entry (insertBackByOuterInner() can not handle empty rows)
      if (columns.size() == 0)
        eigenMatrix_.insertBackByOuterInner(row, row);
    }
    eigenMatrix_.finalize();
    eigenMatrix_.makeCompressed();
  } // RowMajorSparseMatrix(...)

  size_type rows() const
  {
    return eigenMatrix_.rows();
  }

  size_type cols() const
  {
    return eigenMatrix_.cols();
  }

  void add(const size_type i, const size_type j, const ElementType& val)
  {
    eigenMatrix_.coeffRef(i, j) += val;
  }

  void set(const size_type i, const size_type j, const ElementType& val)
  {
    eigenMatrix_.coeffRef(i, j) = val;
  }

  const ElementType get(const size_type i, const size_type j) const
  {
    return eigenMatrix_.coeff(i, j);
  }

  BackendType& backend()
  {
    return eigenMatrix_;
  }

  const BackendType& backend() const
  {
    return eigenMatrix_;
  }

private:
  BackendType eigenMatrix_;
}; // class EigenRowMajorSparseMatrix


template <class ElementImp = double>
class EigenDenseMatrixTraits
{
public:
  typedef ElementImp ElementType;
  typedef EigenDenseMatrix<ElementType> derived_type;
  //  typedef derived_type EigenMatrix_derived_type;
  typedef typename ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, ::Eigen::Dynamic> BackendType;
  typedef typename BackendType::Index size_type;
}; // class DenseMatrixTraits


template <class ElementImp = double>
class EigenDenseMatrix : public MatrixInterface<EigenDenseMatrixTraits<ElementImp>>,
                         public EigenMatrixInterface<EigenDenseMatrixTraits<ElementImp>>
{
public:
  typedef EigenDenseMatrix<ElementImp> ThisType;
  typedef EigenDenseMatrixTraits<ElementImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ElementType ElementType;
  typedef typename Traits::size_type size_type;

  EigenDenseMatrix()
  {
  }

  EigenDenseMatrix(const ThisType& _other)
    : eigenMatrix_(_other.backend())
  {
  }

  EigenDenseMatrix(const BackendType& _otherEigenMatrix)
    : eigenMatrix_(_otherEigenMatrix)
  {
  }

  EigenDenseMatrix(const EigenRowMajorSparseMatrix<ElementType>& _other)
    : eigenMatrix_(_other.backend())
  {
  }

  EigenDenseMatrix(const EigenDenseVector<ElementType>& _other)
    : eigenMatrix_(_other.backend())
  {
  }

  EigenDenseMatrix(const size_type _rows, const size_type _cols)
    : eigenMatrix_(_rows, _cols)
  {
  }

  ThisType& operator=(const ThisType& _other)
  {
    eigenMatrix_ = _other.backend();
    return *this;
  }

  ThisType& operator=(const BackendType& _otherEigenMatrix)
  {
    eigenMatrix_ = _otherEigenMatrix;
    return *this;
  }

  size_type rows() const
  {
    return eigenMatrix_.rows();
  }

  size_type cols() const
  {
    return eigenMatrix_.cols();
  }

  void add(const size_type i, const size_type j, const ElementType& val)
  {
    eigenMatrix_(i, j) += val;
  }

  void set(const size_type i, const size_type j, const ElementType& val)
  {
    eigenMatrix_(i, j) = val;
  }

  const ElementType get(const size_type i, const size_type j) const
  {
    return eigenMatrix_(i, j);
  }

  BackendType& backend()
  {
    return eigenMatrix_;
  }

  const BackendType& backend() const
  {
    return eigenMatrix_;
  }

private:
  BackendType eigenMatrix_;
}; // class EigenDenseMatrix


template <class ElementImp = double>
class EigenDenseVectorTraits
{
public:
  typedef ElementImp ElementType;
  typedef EigenDenseVector<ElementType> derived_type;
  //  typedef derived_type EigenVector_derived_type;
  typedef typename ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, 1> BackendType;
  typedef typename BackendType::Index size_type;
}; // class DenseVectorTraits


template <class ElementImp = double>
class EigenDenseVector : public VectorInterface<EigenDenseVectorTraits<ElementImp>>,
                         public EigenVectorInterface<EigenDenseVectorTraits<ElementImp>>
{
public:
  typedef EigenDenseVector<ElementImp> ThisType;
  typedef EigenDenseVectorTraits<ElementImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ElementType ElementType;
  typedef typename Traits::size_type size_type;

  EigenDenseVector()
  {
  }

  EigenDenseVector(const ThisType& _other)
    : eigenVector_(_other.backend())
  {
  }

  EigenDenseVector(const BackendType& _otherEigenVector)
    : eigenVector_(_otherEigenVector)
  {
  }

  EigenDenseVector(const size_type _size)
    : eigenVector_(_size)
  {
  }

  ThisType& operator=(const ThisType& _other)
  {
    eigenVector_ = _other.backend();
    return *this;
  }

  ThisType& operator=(const BackendType& _otherEigenVector)
  {
    eigenVector_ = _otherEigenVector;
    return *this;
  }

  size_type size() const
  {
    return eigenVector_.size();
  }

  void add(const size_type i, const ElementType& val)
  {
    eigenVector_(i) += val;
  }

  void set(const size_type i, const ElementType& val)
  {
    eigenVector_(i) = val;
  }

  const ElementType get(const size_type i) const
  {
    return eigenVector_(i);
  }

  BackendType& backend()
  {
    return eigenVector_;
  }

  const BackendType& backend() const
  {
    return eigenVector_;
  }

private:
  BackendType eigenVector_;
}; // class DenseVector


template <class ElementType = double, class size_type = unsigned int>
Dune::shared_ptr<EigenRowMajorSparseMatrix<ElementType>> createIdentityEigenRowMajorSparseMatrix(const size_type _size)
{
  // create the sparsity pattern
  SparsityPatternDefault pattern(_size);
  for (typename SparsityPatternDefault::size_type ii = 0; ii < _size; ++ii)
    pattern.set(ii).insert(ii);
  Dune::shared_ptr<EigenRowMajorSparseMatrix<ElementType>> ret =
      Dune::make_shared<EigenRowMajorSparseMatrix<ElementType>>(_size, _size, pattern);
  for (typename SparsityPatternDefault::size_type ii = 0; ii < _size; ++ii)
    ret->set(ii, ii, 1.0);
  return ret;
}


} // namespace Container
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_HH
