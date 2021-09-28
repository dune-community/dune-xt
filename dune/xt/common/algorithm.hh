// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2014 - 2016, 2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_COMMON_ALGORITHM_HH
#define DUNE_XT_COMMON_ALGORITHM_HH

#include <utility>
#include <iterator>
#include <algorithm>

#include <dune/xt/common/string.hh>

namespace Dune::XT::Common {


//! same as std::copy_if, but with move
template <class InputIt, class OutputIt, class UnaryPredicate>
OutputIt move_if(InputIt first, InputIt last, OutputIt d_first, UnaryPredicate pred)
{
  while (first != last) {
    if (pred(*first))
      *d_first++ = std::move(*first);
    first++;
  }
  return d_first;
}

template <class InputIt, class SequenceType = std::vector<std::string>>
SequenceType make_string_sequence(InputIt first, InputIt last)
{
  SequenceType ret;
  const auto to_str = [](const typename std::iterator_traits<InputIt>::value_type& val) { return to_string(val); };
  std::transform(first, last, std::back_inserter(ret), to_str);
  return ret;
}


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_ALGORITHM_HH
