// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2010, 2012 - 2017)
//   René Fritze     (2009 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2011 - 2012)
//   Tobias Leibner  (2014, 2018, 2020)

#ifndef DUNE_XT_COMMON_PRINT_HH
#define DUNE_XT_COMMON_PRINT_HH

#include <string>
#include <iostream>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/filesystem.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Common {
namespace internal {


/// \note Should be used to derive from (except for vectors and matrices), when specializing Printer.
/// \sa Printer
/// \sa VectorPrinter
template <class T, bool use_repr = false>
class DefaultPrinter
{
  using ThisType = DefaultPrinter;

public:
  using ValueType = T;

  const ValueType& value;
  const Configuration opts;

  DefaultPrinter(const ValueType& val, const Configuration& cfg)
    : value(val)
    , opts(cfg)
  {}

  DefaultPrinter(const ThisType&) = default;
  DefaultPrinter(ThisType&&) = default;

  virtual ~DefaultPrinter() = default;

  virtual void repr(std::ostream& out) const
  {
    if constexpr (is_printable<T>::value) {
      out << value;
    } else {
      out << "missing specialization for Printer<T> with T=" << Typename<T>::value();
    }
  }

  virtual void str(std::ostream& out) const
  {
    this->repr(out);
  }

  void print(std::ostream& out) const
  {
    if (use_repr)
      this->repr(out);
    else
      this->str(out);
  }
}; // class DefaultPrinter


/// \note Should be used to derive from when specializing Printer for vectors.
/// \sa Printer
/// \sa DefaultPrinter
template <class T, bool use_repr>
class VectorPrinter : public internal::DefaultPrinter<T, use_repr>
{
public:
  const std::string class_name;

  VectorPrinter(const T& val, const Configuration& cfg = {}, const std::string& clss_nm = Typename<T>::value())
    : internal::DefaultPrinter<T, use_repr>(val, cfg)
    , class_name(clss_nm)
  {}

  void repr(std::ostream& out) const override
  {
    out << class_name << "(";
    const auto sz = this->value.size();
    if (sz > 0) {
      const std::string delim =
          (std::use_facet<std::numpunct<char>>(std::cout.getloc()).decimal_point() == ',') ? ";" : ",";
      out << "{" << print(this->value[0], this->opts) /*V::get_entry(this->value, 0)*/;
      for (auto&& ii : value_range(decltype(sz)(1), sz))
        out << delim << " " << print(this->value[ii], this->opts) /*V::get_entry(this->value, ii)*/;
      out << "}";
    }
    out << ")";
  } // ... repr(...)

  void str(std::ostream& out) const override
  {
    const auto sz = this->value.size();
    if (sz == 0)
      out << "[]";
    else {
      const std::string delim =
          (std::use_facet<std::numpunct<char>>(std::cout.getloc()).decimal_point() == ',') ? ";" : ",";
      out << "[" << print(this->value[0], this->opts) /*V::get_entry(this->value, 0)*/;
      for (auto&& ii : value_range(decltype(sz)(1), sz))
        out << delim << " " << print(this->value[ii], this->opts) /*V::get_entry(this->value, ii)*/;
      out << "]";
    }
  } // ... str(...)
}; // class VectorPrinter


/// \note Should be used to derive from when specializing Printer for matrices.
/// \sa Printer
/// \sa DefaultPrinter
template <class T, bool use_repr>
class MatrixPrinter : public internal::DefaultPrinter<T, use_repr>
{
  static_assert(is_matrix<T>::value, "");

  using M = MatrixAbstraction<T>;

public:
  const std::string class_name;

  MatrixPrinter(const T& val, const Configuration& cfg = {}, const std::string& clss_nm = Typename<T>::value())
    : internal::DefaultPrinter<T, use_repr>(val, cfg)
    , class_name(clss_nm)
  {}

  void repr(std::ostream& out) const override
  {
    out << class_name << "(";
    const auto rows = M::rows(this->value);
    const auto cols = M::cols(this->value);
    if (rows * cols > 0) {
      out << "{";
      const std::string delim =
          (std::use_facet<std::numpunct<char>>(std::cout.getloc()).decimal_point() == ',') ? ";" : ",";
      const std::string newline = "\n";
      for (auto&& ii : value_range(rows)) {
        out << (ii == 0 ? "{" : " ") << "{" << print(M::get_entry(this->value, ii, 0), this->opts);
        for (auto&& jj : value_range(decltype(cols)(1), cols))
          out << delim << " " << print(M::get_entry(this->value, ii, jj), this->opts);
        out << "}" << ((ii == rows - 1) ? "" : ",") << ((ii == rows - 1) ? "" : newline);
      }
      out << "}";
    }
    out << ")";
  } // ... repr(...)

  void str(std::ostream& out) const override
  {
    const auto rows = M::rows(this->value);
    const auto cols = M::cols(this->value);
    out << "[";
    if (rows * cols > 0) {
      const std::string delim =
          (std::use_facet<std::numpunct<char>>(std::cout.getloc()).decimal_point() == ',') ? ";" : ",";
      const std::string newline = this->opts.get("oneline", false) ? "" : "\n";
      for (auto&& ii : value_range(rows)) {
        out << (ii == 0 ? "" : " ") << "[" << print(M::get_entry(this->value, ii, 0), this->opts);
        for (auto&& jj : value_range(decltype(cols)(1), cols))
          out << delim << " " << print(M::get_entry(this->value, ii, jj), this->opts);
        out << "]" << ((ii == rows - 1) ? "" : ",") << ((ii == rows - 1) ? "" : newline);
      }
    }
    out << "]";
  } // ... str(...)
}; // class MatrixPrinter


} // namespace internal


/**
 * \brief Specialize this class for your type to benefit from the print() and repr() functions.
 *
 * Any specialization of Printer is expected to have the following print() method
<code>
  void print(std::ostream& out) const
  {
    if (use_repr)
      this->repr(out);
    else
      this->str(out);
  }
</code>
 * the behaviour of which is supposed to change with use_repr. A possible implemenation is given in
 * internal::DefaultPrinter, based on the following two methods:
<code>
  void repr(std::ostream& out) const
  {
    out << "detailed and unambiguous description";
  }

  void str(std::ostream& out) const
  {
    out << "compact and readabale description";
  }
</code>
 * It is therefore recommended to derive from internal::DefaultPrinter or internal::VectorPrinter, where some of this
 * is already provided.
 *
 * Similar to python, one can use repr(obj) to obtain an unambigous representation of obj (with as much detail as
 * possible), while print(obj) mimics Pythons __str__ to obtain a compact and readable representation of obj.
 *
 * \sa internal::DefaultPrinter
 * \sa print
 * \sa repr
 */
template <class T, bool use_repr = false, typename anything = void>
class Printer : public internal::DefaultPrinter<T, use_repr>
{
public:
  Printer(const T& val, const Configuration& param = {})
    : internal::DefaultPrinter<T, use_repr>(val, param)
  {}
};


/// Specialization of Printer for all our vectors
template <class V, bool use_repr>
class Printer<V, use_repr, std::enable_if_t<is_vector<V>::value>> : public internal::VectorPrinter<V, use_repr>
{
public:
  Printer(const V& val, const Configuration& param)
    : internal::VectorPrinter<V, use_repr>(val, param)
  {}
}; // class Printer


/// Specialization of Printer for all our matrices
template <class M, bool use_repr>
class Printer<M, use_repr, std::enable_if_t<is_matrix<M>::value>> : public internal::MatrixPrinter<M, use_repr>
{
public:
  Printer(const M& val, const Configuration& param = {{"oneline", "false"}})
    : internal::MatrixPrinter<M, use_repr>(val, param)
  {}
}; // class Printer


/// \sa Printer
template <class T, bool use_repr>
std::ostream& operator<<(std::ostream& out, const Printer<T, use_repr>& printer)
{
  printer.print(out);
  return out;
}


/**
 * The purpose of print(obj) is to provide a flexible and customizable means of modifying the behaviour of printing an
 * object to a stream, without conflicting with implementations of operator<< for obj, which are often beyond the reach
 * of the user. Can be used as in:
<code>
std::cout << print(complicated_object, {"tabular", "true"}) << std::endl;
</code>
 * The output can be determined by specializing Printer<TypeOfComplicatedObject>, which keeps a reference to
 * complicated_object and a Configuration({"tabular", "true"}), which can be used by the implementor to customize the
 * output. Defaults to
<code>
std::cout << complicated_object << std::endl;
</code>
 * if operator<< is available for complicated_object and no specialization of Printer is found.
 *
 * \sa Printer
 * \sa repr
 */
template <class T>
Printer<T, false> print(const T& value, const Configuration& param = {})
{
  return Printer<T, false>(value, param);
}


/**
 * The purpose of repr(obj) is to provide a flexible and customizable means of modifying the behaviour of printing an
 * object to a stream, without conflicting with implementations of operator<< for obj, which are often beyond the reach
 * of the user. Can be used as in:
<code>
std::cout << print(complicated_object, {"tabular", "true"}) << std::endl;
</code>
 * The output can be determined by specializing Printer<TypeOfComplicatedObject>, which keeps a reference to
 * complicated_object and a Configuration({"tabular", "true"}), which can be used by the implementor to customize the
 * output. Defaults to
<code>
std::cout << complicated_object << std::endl;
</code>
 * if operator<< is available for complicated_object and no specialization of Printer is found.
 *
 * \sa Printer
 * \sa print
 */
template <class T>
Printer<T, true> repr(const T& value, const Configuration& param = {})
{
  return Printer<T, true>(value, param);
}


//! useful for visualizing sparsity patterns of matrices
template <class Matrix>
void matrix_to_gnuplot_stream(const Matrix& matrix, std::ostream& stream)
{
  unsigned long nz = 0;

  const auto cols = matrix.cols();
  for (auto row : value_range(matrix.rows())) {
    for (auto col : value_range(cols)) {
      if (matrix.find(row, col))
        stream << row << "\t" << col << "\t" << matrix(row, col) << "\n";
    }
    nz += matrix.numNonZeros(row);
    stream << "#non zeros in row " << row << " " << matrix.numNonZeros(row) << " (of " << matrix.cols() << " cols)\n";
  }
  stream << "#total non zeros " << nz << " of " << matrix.rows() * matrix.cols() << " entries\n";
} // matrix_to_gnuplot_stream


//! maps 1,2,3 to x,y,z / X,Y,Z
inline std::string dim_to_axis_name(const size_t dim, const bool capitalize = false)
{
  char c = 'x';

  c += dim;
  if (capitalize)
    c -= 32;
  return std::string() += c;
} // matrix_to_gnuplot_stream


/** Outputiterator to emulate python's str.join(iterable)
 * \see http://codereview.stackexchange.com/questions/30132/comma-formatted-stl-vectors/30181#30181
 * \example std::copy(strings.begin(), strings.end(), PrefixOutputIterator<string>(ostream, ","));
 **/
template <typename T>
class PrefixOutputIterator
{
  std::ostream& ostream;
  std::string prefix;
  bool first;

public:
  using difference_type = std::size_t;
  using value_type = T;
  using pointer = T*;
  using reference = T;
  using iterator_category = std::output_iterator_tag;

  PrefixOutputIterator(std::ostream& o, std::string const& p = "")
    : ostream(o)
    , prefix(p)
    , first(true)
  {}

  PrefixOutputIterator& operator*()
  {
    return *this;
  }
  PrefixOutputIterator& operator++()
  {
    return *this;
  }
  PrefixOutputIterator& operator++(int)
  {
    return *this;
  }

  void operator=(T const& value)
  {
    if (first) {
      ostream << value;
      first = false;
    } else {
      ostream << prefix << value;
    }
  }
}; // class PrefixOutputIterator


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_PRINT_HH
