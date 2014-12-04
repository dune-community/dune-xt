// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH
#define DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH

#ifndef DUNE_STUFF_SSIZE_T
#define DUNE_STUFF_SSIZE_T long int
#endif

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
 *  \brief  Contains tags mostly needed for python bindings.
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

  /**
   * \defgroup haveto ´´These methods have to be implemented by a derived class in addition to the ones required by
   * ContainerInterface!``
   * \{
   */

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

  template <class XX>
  typename XX::derived_type operator*(const VectorInterface<XX, ScalarType>& xx) const
  {
    typename XX::derived_type yy(cols());
    mv(xx.as_imp(xx), yy);
    return yy;
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

  /**
   * \}
   */

  ScalarType sup_norm() const
  {
    ScalarType ret = 0;
    for (size_t ii = 0; ii < rows(); ++ii)
      for (size_t jj = 0; jj < cols(); ++jj)
        ret = std::max(ret, std::abs(get_entry(ii, jj)));
    return ret;
  }

  /**
   * \defgroup python_bindings ´´These methods are necesarry for the python bindings. They are provided by the
   * interface!``
   * \{
   */

  inline DUNE_STUFF_SSIZE_T pb_rows() const
  {
    if (!(rows() <= std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()))
      DUNE_THROW(Exceptions::index_out_of_range,
                 "The number of rows of this (" << rows() << ") do not fit into DUNE_STUFF_SSIZE_T (max "
                                                << std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()
                                                << ")!");
    return (DUNE_STUFF_SSIZE_T)(rows());
  } // ... pb_rows(...)

  inline DUNE_STUFF_SSIZE_T pb_cols() const
  {
    if (!(cols() <= std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()))
      DUNE_THROW(Exceptions::index_out_of_range,
                 "The number of columns of this (" << cols() << ") do not fit into DUNE_STUFF_SSIZE_T (max "
                                                   << std::numeric_limits<DUNE_STUFF_SSIZE_T>::max()
                                                   << ")!");
    return (DUNE_STUFF_SSIZE_T)(cols());
  } // ... pb_cols(...)

  inline void pb_add_to_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    add_to_entry(
        this->assert_is_size_t_compatible_and_convert(ii), this->assert_is_size_t_compatible_and_convert(jj), value);
  } // ... pb_add_to_entry(...)

  inline void pb_set_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj, const ScalarType& value)
  {
    set_entry(
        this->assert_is_size_t_compatible_and_convert(ii), this->assert_is_size_t_compatible_and_convert(jj), value);
  } // ... pb_set_entry(...)

  inline ScalarType pb_get_entry(const DUNE_STUFF_SSIZE_T ii, const DUNE_STUFF_SSIZE_T jj) const
  {
    return get_entry(this->assert_is_size_t_compatible_and_convert(ii),
                     this->assert_is_size_t_compatible_and_convert(jj));
  } // ... pb_get_entry(...)

  inline void pb_clear_row(const size_t ii)
  {
    clear_row(this->assert_is_size_t_compatible_and_convert(ii));
  }

  inline void pb_clear_col(const size_t jj)
  {
    clear_col(this->assert_is_size_t_compatible_and_convert(jj));
  }

  inline void pb_unit_row(const size_t ii)
  {
    unit_row(this->assert_is_size_t_compatible_and_convert(ii));
  }

  inline void pb_unit_col(const size_t jj)
  {
    unit_col(this->assert_is_size_t_compatible_and_convert(jj));
  }

  /**
   * \}
   */

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


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_CONTAINER_MATRIX_INTERFACE_HH
