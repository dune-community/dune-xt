#ifndef DUNE_STUFF_LA_CONTAINER_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_INTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune {
namespace Stuff {
namespace LA {


template <class Traits>
class ContainerInterface
{
public:
  typedef typename Traits::BackendType BackendType;

  virtual BackendType& backend() = 0;
  virtual ~ContainerInterface()
  {
  }

  virtual const BackendType& backend() const = 0;
}; // class Interface


template <class Traits>
class MatrixInterface : public ContainerInterface<Traits>
{
public:
  typedef MatrixInterface<Traits> ThisType;
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ElementType ElementType;

  size_t rows() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().rows());
    return asImp().rows();
  }

  size_t cols() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().cols());
    return asImp().cols();
  }

  void add(const size_t i, const size_t j, const ElementType& val)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().add(i, j, val));
  }

  void set(const size_t i, const size_t j, const ElementType& val)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().set(i, j, val));
  }

  const ElementType get(const size_t i, const size_t j) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().get(i, j));
    return asImp().get(i, j);
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class MatrixInterface


template <class Traits>
class VectorInterface : public ContainerInterface<Traits>
{
public:
  typedef VectorInterface<Traits> ThisType;
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::ElementType ElementType;

  size_t size() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  void add(const size_t i, const ElementType& val)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().add(i, val));
  }

  void set(const size_t i, const ElementType& val)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().set(i, val));
  }

  const ElementType get(const size_t i) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().get(i));
    return asImp().get(i);
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class VectorInterface


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_INTERFACE_HH
