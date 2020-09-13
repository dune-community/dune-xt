// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2009, 2012 - 2017)
//   René Fritze     (2009 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2011 - 2012)
//   Tobias Leibner  (2014, 2018, 2020)

#ifndef DUNE_XT_COMMON_MISC_HH
#define DUNE_XT_COMMON_MISC_HH

#include <cstring>
#include <map>
#include <assert.h>
#include <algorithm>
#include <dune/common/version.hh>
#include <dune/common/exceptions.hh>

#include <cstddef>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cerrno>
#include <limits>
#include <string.h>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include <dune/xt/common/disable_warnings.hh>
#include <boost/algorithm/string.hpp>
#include <dune/xt/common/reenable_warnings.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>

namespace Dune {
namespace XT {
namespace Common {


//! element-index-in-container search
template <class StlSequence>
inline int get_idx(const StlSequence& ct, const typename StlSequence::value_type& val)
{
  const auto result = std::find(ct.begin(), ct.end(), val);
  if (result == ct.end())
    return -1;
  return std::distance(ct.begin(), result);
} // get_idx

//! type safe (this will not compile for degraded-to-pointer arrays) way of getting array length
template <class T, size_t N>
size_t array_length(T (&/*array*/)[N])
{
  return N;
}

//! get a non-zero initialised std::array
template <class T, size_t N>
std::array<T, N> make_array(const T& v)
{
  std::array<T, N> ret;
  ret.fill(v);
  return ret;
}

template <class T, size_t N>
std::array<T, N> make_array(const std::vector<T>& v)
{
  if (v.size() == 1) {
    return make_array<T, N>(v[0]);
  } else if (v.size() < N) {
    std::array<T, N> error_msg;
    DUNE_THROW(Exceptions::shapes_do_not_match,
               "Trying to create a " << get_typename(error_msg) << " from a " << get_typename(v) << " of size "
                                     << v.size());
  } else {
    std::array<T, N> ret;
    for (size_t nn = 0; nn < N; ++nn) {
      ret[nn] = v[nn];
    }
    return ret;
  }
}

//! writes process environment to file
void dump_environment(boost::filesystem::ofstream& file, std::string csv_sep = ",");

} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_MISC_HH
