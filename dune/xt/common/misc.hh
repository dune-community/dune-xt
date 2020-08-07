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
#include <boost/filesystem/fstream.hpp>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/deprecated.hh>

namespace Dune::XT::Common {

//! element-index-in-container search
template <class StlSequence>
inline int get_idx(const StlSequence& ct, const typename StlSequence::value_type& val)
{
  const auto result = std::find(ct.begin(), ct.end(), val);
  if (result == ct.end())
    return -1;
  return std::distance(ct.begin(), result);
} // get_idx

/** this allows subscription indices to wrap around
 * \example N=4: WraparoundArray[4] == WraparoundArray[0] && WraparoundArray[-1] == WraparoundArray[3]
 **/
template <class T, size_t N>
struct WraparoundArray : public std::array<T, N>
{
  typedef std::array<T, N> BaseType;
  WraparoundArray()
  {
    for (size_t i = 0; i < N; ++i)
      this->operator[](i) = T();
  }

  WraparoundArray(const BaseType other)
  {
    for (size_t i = 0; i < N; ++i)
      this->operator[](i) = other[i];
  }

  typename BaseType::reference operator[](std::size_t i)
  {
    return BaseType::operator[](i % N);
  }

  typename BaseType::reference operator[](int i)
  {
    const std::size_t real_index = i < 0 ? static_cast<std::size_t>(N - (((i * -1) % N) + 1)) : std::size_t(i);
    return BaseType::operator[](real_index);
  } // []

  typename BaseType::const_reference operator[](std::size_t i) const
  {
    return BaseType::operator[](i % N);
  }

  typename BaseType::const_reference operator[](int i) const
  {
    const std::size_t real_index = i < 0 ? static_cast<std::size_t>(N - (((i * -1) % N) + 1)) : std::size_t(i);
    return BaseType::operator[](real_index);
  } // []
};

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

/** abstraction for stdlibs w/o map.emplace
 * This isn't covering the full emplace spec, and the piecewise case potentially improper,
 * but we only need this in an Intel MIC setup with weird lib versioning.
 * Normally we'd depend on the full gcc 4.8 stack atm and do not use the insert fallback
 **/
template <typename Key, typename T, typename MapType>
std::pair<typename MapType::iterator, bool> map_emplace(MapType& map_in, Key key, T value)
{
#if HAVE_MAP_EMPLACE
  return map_in.emplace(key, value);
#else
  DXTC_LOG_DEBUG_0 << "using map.insert fallback instead of emplace\n";
  return map_in.insert(typename MapType::value_type(key, value));
#endif
}

template <typename K, typename V, typename MapType>
DXT_DEPRECATED_MSG("no longer needed with c++14 (2018/03/20)")
std::pair<typename MapType::iterator, bool> map_emplace(MapType& map_in,
                                                        std::piecewise_construct_t pcw,
                                                        K&& keys,
                                                        V&& values)
{
  return map_in.emplace(pcw, keys, values);
}

template <typename T>
struct remove_const_reference
{
  using type = typename std::remove_reference<typename std::remove_const<T>::type>::type;
};

} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_MISC_HH
