// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#ifndef PYTHON_DUNE_XT_COMMON_NUMPY_HH
#define PYTHON_DUNE_XT_COMMON_NUMPY_HH

#include <dune/pybindxi/numpy.h>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/numeric.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/string.hh>

namespace Dune::XT::Common::bindings {


template <class T>
std::string print_shape(const pybind11::array_t<T>& array)
{
  std::string ret = "(";
  if (array.ndim() > 0) {
    ret += to_string(array.shape(0));
    for (ssize_t nn = 1; nn < array.ndim(); ++nn)
      ret += ", " + to_string(array.shape(nn));
  }
  ret += ")";
  return ret;
}


template <size_t ndim, class T>
auto access_array(const pybind11::array_t<T>& array,
                  const std::array<size_t, ndim> required_shape,
                  const std::string& array_name = "",
                  const std::string& required_shape_docs = "")
{
  // check
  bool all_good = true;
  if (numeric_cast<size_t>(array.ndim()) != ndim)
    all_good = false;
  for (size_t ii = 0; ii < ndim; ++ii)
    if (numeric_cast<size_t>(array.shape(ii)) != required_shape[ii])
      all_good = false;
  // report
  DUNE_THROW_IF(!all_good,
                XT::Common::Exceptions::shapes_do_not_match,
                (array_name.empty() ? "numpy.ndarray" : array_name)
                    << " has wrong shape!"
                    << "\n   "
                    << "should be "
                    << (required_shape_docs.empty() ? to_string(required_shape)
                                                    : required_shape_docs + ": " + to_string(required_shape))
                    << "\n   "
                    << "is " << print_shape(array) << "!");
  // return access
  return array.template unchecked<ndim>();
} // ... access_array(...)


template <size_t ndim, class T>
auto access_array(pybind11::array_t<T>& array,
                  const std::array<size_t, ndim> required_shape,
                  const std::string& array_name = "",
                  const std::string& required_shape_docs = "")
{
  // check
  bool all_good = true;
  if (numeric_cast<size_t>(array.ndim()) != ndim)
    all_good = false;
  for (size_t ii = 0; ii < ndim; ++ii)
    if (numeric_cast<size_t>(array.shape(ii)) != required_shape[ii])
      all_good = false;
  // report
  DUNE_THROW_IF(!all_good,
                XT::Common::Exceptions::shapes_do_not_match,
                (array_name.empty() ? "numpy.ndarray" : array_name)
                    << " has wrong shape!"
                    << "\n   "
                    << "should be "
                    << (required_shape_docs.empty() ? to_string(required_shape)
                                                    : required_shape_docs + ": " + to_string(required_shape))
                    << "\n   "
                    << "is " << print_shape(array) << "!");
  // return access
  return array.template mutable_unchecked<ndim>();
} // ... access_array(...)


} // namespace Dune::XT::Common::bindings

#endif // PYTHON_DUNE_XT_COMMON_NUMPY_HH
