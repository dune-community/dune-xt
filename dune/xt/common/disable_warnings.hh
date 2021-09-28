// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2014, 2016 - 2017)
//   René Fritze     (2011 - 2012, 2014 - 2019)
//   Tobias Leibner  (2014, 2016 - 2020)

#include <boost/config.hpp>

// Please sort these alphabetically!
#if defined(BOOST_CLANG) && BOOST_CLANG
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wcovered-switch-default"
#  pragma clang diagnostic ignored "-Wdelete-non-virtual-dtor"
#  pragma GCC diagnostic ignored "-Wdeprecated"
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#  pragma GCC diagnostic ignored "-Wdeprecated-register"
#  pragma GCC diagnostic ignored "-Wdisabled-macro-expansion"
#  pragma GCC diagnostic ignored "-Wextra"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#  pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#  pragma GCC diagnostic ignored "-Wignored-qualifiers"
#  pragma GCC diagnostic ignored "-Wlogical-not-parentheses"
#  pragma GCC diagnostic ignored "-Wlogical-op-parentheses"
#  pragma GCC diagnostic ignored "-Wmismatched-tags"
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Woverloaded-virtual"
#  pragma GCC diagnostic ignored "-Wpedantic"
#  pragma GCC diagnostic ignored "-Wredundant-decls"
#  pragma GCC diagnostic ignored "-Wreorder"
#  pragma GCC diagnostic ignored "-Wreturn-stack-address"
#  pragma GCC diagnostic ignored "-Wshadow"
#  pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#  pragma GCC diagnostic ignored "-Wswitch-enum"
#  pragma GCC diagnostic ignored "-Wtype-limits"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wuninitialized"
#  pragma GCC diagnostic ignored "-Wunreachable-code"
#  pragma GCC diagnostic ignored "-Wunused-const-variable"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wunused-private-field"
#  pragma GCC diagnostic ignored "-Wused-but-marked-unused"
#elif defined(BOOST_GCC) && BOOST_GCC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wattributes"
#  if (__GNUC__ > 8)
#    pragma GCC diagnostic ignored "-Wdeprecated-copy"
#  endif
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#  pragma GCC diagnostic ignored "-Wextra"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#  pragma GCC diagnostic ignored "-Wignored-qualifiers"
#  pragma GCC diagnostic ignored "-Wlogical-op"
#  if (__GNUC__ > 5)
#    pragma GCC diagnostic ignored "-Wmisleading-indentation"
#  endif
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Woverloaded-virtual"
#  pragma GCC diagnostic ignored "-Wpedantic"
#  pragma GCC diagnostic ignored "-Wredundant-decls"
#  pragma GCC diagnostic ignored "-Wreorder"
#  pragma GCC diagnostic ignored "-Wshadow"
#  pragma GCC diagnostic ignored "-Wswitch-enum"
#  pragma GCC diagnostic ignored "-Wtype-limits"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wuninitialized"
#  pragma GCC diagnostic ignored "-Wunreachable-code"
#  pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#elif defined(BOOST_INTEL) && BOOST_INTEL
#  pragma warning push
#  pragma warning disable 177 187 654 1011 2259 2282
#endif
