#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#if HAVE_EIGEN

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <dune/common/bartonnackmanifcheck.hh>

#include "interface.hh"
#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {
namespace Eigen {


template <class Traits>
class Base
{
public:
  typedef typename Traits::derived_type derived_type;
}; // class Base


template <class Traits>
class MatrixBase : public Base<Traits>
{
public:
  typedef typename Base<Traits>::derived_type derived_type;
}; // class MatrixBase


template <class Traits>
class VectorBase : public Base<Traits>
{
public:
  typedef typename Base<Traits>::derived_type derived_type;
}; // class VectorBase


template <class ElementType>
class RowMajorSparseMatrix;


template <class ElementImp = double>
class RowMajorSparseMatrixTraits
{
public:
  typedef ElementImp ElementType;

  typedef RowMajorSparseMatrix<ElementType> derived_type;

  typedef typename ::Eigen::SparseMatrix<ElementType, ::Eigen::RowMajor> BaseType;

  typedef typename BaseType::Index size_type;
}; // class RowMajorSparseMatrixTraits


template <class ElementImp = double>
class RowMajorSparseMatrix : public MatrixInterface<RowMajorSparseMatrixTraits<ElementImp>>,
                             virtual public RowMajorSparseMatrixTraits<ElementImp>::BaseType,
                             public MatrixBase<RowMajorSparseMatrixTraits<ElementImp>>
{
public:
  typedef RowMajorSparseMatrixTraits<ElementImp> Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

  RowMajorSparseMatrix()
    : BaseType()
  {
  }

  RowMajorSparseMatrix(const size_type _rows, const size_type _cols,
                       const Dune::Stuff::LA::Container::Pattern::Default& _pattern)
    : BaseType(_rows, _cols)
  {
    assert(size_type(_pattern.rows()) == _rows && "Given pattern too short!");
    typedef Dune::Stuff::LA::Container::Pattern::Default PatternType;
    typedef PatternType::ColumnsType ColumnsType;
    for (size_type row = 0; row < size_type(_pattern.rows()); ++row) {
      BaseType::startVec(row);
      const ColumnsType& columns = _pattern.columns(row);
      for (typename ColumnsType::const_iterator columnIt = columns.begin(); columnIt != columns.end(); ++columnIt) {
        const size_type column = *columnIt;
        BaseType::insertBackByOuterInner(row, column);
      }
      // create diagonal entry (insertBackByOuterInner() can not handle empty rows)
      if (columns.size() == 0)
        BaseType::insertBackByOuterInner(row, row);
    }
    BaseType::finalize();
    BaseType::makeCompressed();
  } // RowMajorSparseMatrix(...)

  size_type rows() const
  {
    return BaseType::rows();
  }

  size_type cols() const
  {
    return BaseType::cols();
  }

  void add(const size_type i, const size_type j, const ElementType& val)
  {
    BaseType::coeffRef(i, j) += val;
  }

  void set(const size_type i, const size_type j, const ElementType& val)
  {
    BaseType::coeffRef(i, j) = val;
  }

  const ElementType get(const size_type i, const size_type j) const
  {
    return BaseType::coeff(i, j);
  }

  BaseType& base()
  {
    return *this;
  }

  const BaseType& base() const
  {
    return *this;
  }
}; // class RowMajorSparseMatrix


template <class ElementImp>
class DenseMatrix;


template <class ElementImp = double>
class DenseMatrixTraits
{
public:
  typedef ElementImp ElementType;

  typedef DenseMatrix<ElementType> derived_type;

  typedef typename ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, ::Eigen::Dynamic> BaseType;

  typedef typename BaseType::Index size_type;
}; // class DenseMatrixTraits


template <class ElementImp = double>
class DenseMatrix : public MatrixInterface<DenseMatrixTraits<ElementImp>>,
                    virtual public DenseMatrixTraits<ElementImp>::BaseType,
                    public MatrixBase<DenseMatrixTraits<ElementImp>>
{
public:
  typedef DenseMatrixTraits<ElementImp> Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

  DenseMatrix(const size_type _rows, const size_type _cols)
    : BaseType(_rows, _cols)
  {
  }

  size_type rows() const
  {
    return BaseType::rows();
  }

  size_type cols() const
  {
    return BaseType::cols();
  }

  void add(const size_type i, const size_type j, const ElementType& val)
  {
    BaseType::operator()(i, j) += val;
  }

  void set(const size_type i, const size_type j, const ElementType& val)
  {
    BaseType::operator()(i, j) = val;
  }

  const ElementType get(const size_type i, const size_type j) const
  {
    return BaseType::operator()(i, j);
  }

  BaseType& base()
  {
    return *this;
  }

  const BaseType& base() const
  {
    return *this;
  }
}; // class DenseMatrix


template <class ElementImp>
class DenseVector;


template <class ElementImp = double>
class DenseVectorTraits
{
public:
  typedef ElementImp ElementType;

  typedef DenseVector<ElementType> derived_type;

  typedef typename ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, 1> BaseType;

  typedef typename BaseType::Index size_type;
}; // class DenseVectorTraits


template <class ElementImp = double>
class DenseVector : public VectorInterface<DenseVectorTraits<ElementImp>>,
                    virtual public DenseVectorTraits<ElementImp>::BaseType,
                    public VectorBase<DenseVectorTraits<ElementImp>>
{
public:
  typedef DenseVectorTraits<ElementImp> Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

  DenseVector(const size_type _size)
    : BaseType(_size)
  {
    BaseType::setZero();
  }

  size_type size() const
  {
    return BaseType::size();
  }

  void add(const size_type i, const ElementType& val)
  {
    BaseType::operator()(i) += val;
  }

  void set(const size_type i, const ElementType& val)
  {
    BaseType::operator()(i) = val;
  }

  const ElementType get(const size_type i) const
  {
    return BaseType::operator()(i);
  }

  BaseType& base()
  {
    return *this;
  }

  const BaseType& base() const
  {
    return *this;
  }
}; // class DenseVector

} // namespace Container
} // namespace Eigen
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_HH
