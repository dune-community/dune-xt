#ifndef DUNE_STUFF_LA_CONTAINER_EIGEN_HH
#define DUNE_STUFF_LA_CONTAINER_EIGEN_HH

#if HAVE_EIGEN

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <dune/common/bartonnackmanifcheck.hh>

#include "pattern.hh"

namespace Dune {
namespace Stuff {
namespace LA {
namespace Container {
namespace Eigen {

template <class MatrixImpTraits>
class MatrixInterface : virtual public MatrixImpTraits::BaseType
{
public:
  typedef MatrixInterface<MatrixImpTraits> ThisType;

  typedef MatrixImpTraits Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

  size_type rows() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().rows());
    return asImp().rows();
  }

  size_type cols() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().cols());
    return asImp().cols();
  }

  void add(const size_type i, const size_type j, const ElementType& val)
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().add(i, j, val));
    asImp().add(i, j, val);
  }

  void set(const size_type i, const size_type j, const ElementType& val)
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().set(i, j, val));
    asImp().set(i, j, val);
  }

  const ElementType get(const size_type i, const size_type j) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().get(i, j));
    return asImp().get(i, j);
  }

  BaseType& base()
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().base());
    return asImp().base();
  }

  const BaseType& base() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().base());
    return asImp().base();
  }

private:
  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class MatrixInterface


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
                             virtual public RowMajorSparseMatrixTraits<ElementImp>::BaseType
{
public:
  typedef RowMajorSparseMatrixTraits<ElementImp> Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

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
                    virtual public DenseMatrixTraits<ElementImp>::BaseType
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


template <class VectorImpTraits>
class VectorInterface : virtual public VectorImpTraits::BaseType
{
public:
  typedef VectorInterface<VectorImpTraits> ThisType;

  typedef VectorImpTraits Traits;

  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BaseType BaseType;

  typedef typename Traits::ElementType ElementType;

  typedef typename Traits::size_type size_type;

  size_type size() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().rows());
    return asImp().size();
  }

  void add(const size_type i, const ElementType& val)
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().add(i, val));
    asImp().add(i, val);
  }

  void set(const size_type i, const ElementType& val)
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().set(i, val));
    asImp().set(i, val);
  }

  const ElementType get(const size_type i) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().get(i));
    return asImp().get(i);
  }

  BaseType& base()
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().base());
    return asImp().base();
  }

  const BaseType& base() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().base());
    return asImp().base();
  }

private:
  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class VectorInterface


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
                    virtual public DenseVectorTraits<ElementImp>::BaseType
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
