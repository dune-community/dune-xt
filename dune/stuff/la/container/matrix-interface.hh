// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH

#include <cmath>
#include <limits>
#include <iostream>
#include <type_traits>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>

#include "container-interface.hh"
#include "vector-interface.hh"

namespace Dune {
namespace Stuff {
namespace LA {


/**
 * \brief Contains tags mostly needed for python bindings.
 */
namespace Tags {


class MatrixInterface
{
};


} // namespace Tags


template <class Traits, class ScalarImp = typename Traits::ScalarType>
class MatrixInterface : public ContainerInterface<Traits, ScalarImp>, public Tags::MatrixInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef ScalarImp ScalarType;

  virtual ~MatrixInterface()
  {
  }

  /// \name Have to be implemented by a derived class in addition to the ones required by ContainerInterface!
  /// \{

  inline size_t rows() const
  {
    CHECK_CRTP(this->as_imp().rows());
    return this->as_imp().rows();
  }

  inline size_t cols() const
  {
    CHECK_CRTP(this->as_imp().cols());
    return this->as_imp().cols();
  }

  template <class XX, class YY>
  inline void mv(const XX& xx, YY& yy) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().mv(xx, yy));
  }

  inline void add_to_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().add_to_entry(ii, jj, value));
  }

  inline void set_entry(const size_t ii, const size_t jj, const ScalarType& value)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().set_entry(ii, jj, value));
  }

  inline ScalarType get_entry(const size_t ii, const size_t jj) const
  {
    CHECK_CRTP(this->as_imp().get_entry(ii, jj));
    return this->as_imp().get_entry(ii, jj);
  }

  inline void clear_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().clear_row(ii));
  }

  inline void clear_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().clear_col(jj));
  }

  inline void unit_row(const size_t ii)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().unit_row(ii));
  }

  inline void unit_col(const size_t jj)
  {
    CHECK_AND_CALL_CRTP(this->as_imp().unit_col(jj));
  }

  /**
   * \brief  Checks entries for inf or nan.
   * \return false if any entry is inf or nan, else true
   */
  inline bool valid() const
  {
    CHECK_CRTP(this->as_imp().valid());
    return this->as_imp().valid();
  }

  /// \}
  /// \name Provided by the interface for convenience.
  /// \note Those marked with vitual should be overriden by any devired class that can do better.
  /// \{

  template <class XX>
  typename XX::derived_type operator*(const VectorInterface<XX, ScalarType>& xx) const
  {
    typename XX::derived_type yy(cols());
    mv(xx.as_imp(xx), yy);
    return yy;
  }

  virtual ScalarType sup_norm() const
  {
    ScalarType ret = 0;
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        ret = std::max(ret, std::abs(get_entry(ii, jj)));
    return ret;
  } // ... sup_norm(...)

  /// \}
  /// \name Necesarry for the python bindings.
  /// \{

  inline DUNE_STUFF_SSIZE_T pb_rows() const
  {
    try {
      return boost::numeric_cast<DUNE_STUFF_SSIZE_T>(rows());
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << rows() << "' to '"
                                                            << Common::Typename<ScalarType>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_rows(...)

  inline DUNE_STUFF_SSIZE_T pb_cols() const
  {
    try {
      return boost::numeric_cast<DUNE_STUFF_SSIZE_T>(cols());
    } catch (boost::bad_numeric_cast& ee) {
      DUNE_THROW(Exceptions::external_error,
                 "There was an error in boost converting '" << cols() << "' to '"
                                                            << Common::Typename<ScalarType>::value()
                                                            << "': "
                                                            << ee.what());
    }
  } // ... pb_cols(...)

  inline void pb_add_to_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    add_to_entry(boost::numeric_cast<size_t>(ii), boost::numeric_cast<size_t>(jj), value);
  }

  inline void pb_set_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    set_entry(boost::numeric_cast<size_t>(ii), boost::numeric_cast<size_t>(jj), value);
  }

  inline ScalarType pb_get_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj) const
  {
    return get_entry(boost::numeric_cast<size_t>(ii), boost::numeric_cast<size_t>(jj));
  }

  inline void pb_clear_row(const DUNE_STUFF_SSIZE_T ii)
  {
    clear_row(boost::numeric_cast<size_t>(ii));
  }

  inline void pb_clear_col(const DUNE_STUFF_SSIZE_T jj)
  {
    clear_col(boost::numeric_cast<size_t>(jj));
  }

  inline void pb_unit_row(const DUNE_STUFF_SSIZE_T ii)
  {
    unit_row(boost::numeric_cast<size_t>(ii));
  }

  inline void pb_unit_col(const DUNE_STUFF_SSIZE_T jj)
  {
    unit_col(boost::numeric_cast<size_t>(jj));
  }

  /// \}

private:
  template <class T, class S>
  friend std::ostream& operator<<(std::ostream& /*out*/, const MatrixInterface<T, S>& /*matrix*/);
}; // class MatrixInterface


template <class T, class S>
std::ostream& operator<<(std::ostream& out, const MatrixInterface<T, S>& matrix)
{
  out << "[";
  const size_t rows = matrix.rows();
  const size_t cols = matrix.cols();
  if (rows > 0 && cols > 0) {
    for (size_t ii = 0; ii < rows; ++ii) {
      if (ii > 0)
        out << "\n ";
      out << "[" << matrix.get_entry(ii, 0);
      for (size_t jj = 1; jj < cols; ++jj)
        out << " " << matrix.get_entry(ii, jj);
      out << "]";
      if (rows > 1 && ii < (rows - 1))
        out << ",";
    }
    out << "]";
  } else
    out << "[ ]]";
  return out;
} // ... operator<<(...)


namespace internal {


template <class M>
struct is_matrix_helper
{
  DSC_has_typedef_initialize_once(Traits) DSC_has_typedef_initialize_once(ScalarType)

      static const bool is_candidate = DSC_has_typedef(Traits)<M>::value && DSC_has_typedef(ScalarType)<M>::value;
}; // class is_matrix_helper


} // namespace internal


template <class M, bool candidate = internal::is_matrix_helper<M>::is_candidate>
struct is_matrix : public std::is_base_of<MatrixInterface<typename M::Traits, typename M::ScalarType>, M>
{
};


template <class M>
struct is_matrix<M, false> : public std::false_type
{
};


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH
