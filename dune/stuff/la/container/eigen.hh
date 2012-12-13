#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#if HAVE_EIGEN

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {
namespace Eigen {

template <class ElementImp = double>
class SparseMatrix : ::Eigen::SparseMatrix<ElementImp, ::Eigen::RowMajor>
{
public:
  typedef ElementImp ElementType;

  typedef ::Eigen::SparseMatrix<ElementType, ::Eigen::RowMajor> BaseType;

  typedef SparseMatrix<ElementType> ThisType;

  typedef typename BaseType::Index size_type;

  SparseMatrix(const size_type _rows, const size_type _cols)
    : BaseType(_rows, _cols)
  {
  }

  SparseMatrix(const size_type _rows, const size_type _cols,
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
    BaseType::setZero();
  } // SparseMatrix(...)

  SparseMatrix(const size_type _rows, const size_type _cols, const size_type _nonZerosPerRow)
    : BaseType(_rows, _cols)
  {
    BaseType::reserve(_nonZerosPerRow);
    BaseType::setZero();
  }

  void reserve(const size_type _nonZerosPerRow)
  {
    BaseType::reserve(_nonZerosPerRow);
    BaseType::setZero();
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
}; // class SparseMatrix


template <class ElementImp = double>
class DenseMatrix : ::Eigen::Matrix<ElementImp, ::Eigen::Dynamic, ::Eigen::Dynamic>
{
public:
  typedef ElementImp ElementType;

  typedef ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, ::Eigen::Dynamic> BaseType;

  typedef DenseMatrix<ElementType> ThisType;

  typedef typename BaseType::Index size_type;

  DenseMatrix(const size_type _rows, const size_type _cols)
    : BaseType(_rows, _cols)
  {
    BaseType::setZero();
  }

  DenseMatrix(const size_type _rows, const size_type _cols, const Dune::Stuff::LA::Container::Pattern::Default&)
    : BaseType(_rows, _cols)
  {
    BaseType::setZero();
  } // SparseMatrix(...)

  DenseMatrix(const size_type _rows, const size_type _cols, const size_type)
    : BaseType(_rows, _cols)
  {
    BaseType::setZero();
  }

  void reserve(const size_type)
  {
    BaseType::setZero();
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


template <class ElementImp = double>
class DenseVector : public ::Eigen::Matrix<ElementImp, ::Eigen::Dynamic, 1>
{
public:
  typedef ElementImp ElementType;

  typedef typename ::Eigen::Matrix<ElementType, ::Eigen::Dynamic, 1> BaseType;

  typedef typename BaseType::Index size_type;

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

  const ElementType& operator[](const size_type i) const
  {
    return BaseType::coeff(i);
  }

  ElementType& operator[](const size_type i)
  {
    return BaseType::coeffRef(i);
  }
}; // class DenseVector

} // namespace Container
} // namespace Eigen
} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_STUFF_LA_CONTAINER_EIGEN_HH
