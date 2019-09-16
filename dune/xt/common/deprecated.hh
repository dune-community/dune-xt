// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze (2018 - 2019)

#ifndef DUNE_XT_COMMON_DEPRECATED_HH
#define DUNE_XT_COMMON_DEPRECATED_HH

#include <dune/common/deprecated.hh>

/** usage examples

 // Deprecate a function
DXT_DEPRECATED
void foo();

// Deprecate a variable
DXT_DEPRECATED
int x;

// Deprecate one declarator in a multi-declarator declaration
int y DXT_DEPRECATED, z;

// Deprecate a function parameter
int triple(DXT_DEPRECATED int x);

// Deprecate a class (or struct)
class DXT_DEPRECATED my_class {
  public:
    // Deprecate a member
    DXT_DEPRECATED int member;
};

// Deprecate an enum
enum DXT_DEPRECATED animals {
  CAT, DOG, MOUSE
};

// Deprecate a typedef
DXT_DEPRECATED
typedef int type;

// Deprecate a template specialization
template <typename T> class templ;

template <>
class DXT_DEPRECATED templ<int> {};

 */

#define DXT_DEPRECATED [[deprecated]]
#define DXT_DEPRECATED_MSG(text) [[deprecated(#text)]]

#endif // DUNE_XT_COMMON_DEPRECATED_HH
