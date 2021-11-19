```
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013 - 2014, 2016 - 2017, 2019)
#   René Fritze     (2009, 2015, 2016, 2019)
#   Tobias Leibner  (2021)
```

# Help

If you are experiencing problems or would like to get help, just send us an [email](dune-xt@dune-community.ovh).
Please note that dune-xt is not one of the core DUNE modules. Thus, you will not get
any support from the official channels (such as the DUNE mailinglist, the DUNE bugtracker,
etc...).
When submitting bugs, please read these
[general guidelines](https://www.dune-project.org/doc/guides/bug_reporting/) beforehand.

# Contributing

We very much encourage and welcome contributions!
There are several ways to contribute to this project:

## Contribute code and/or documentation

Fork the repository and open up a pull request, explaining your changes and why they should be included.
Before you start making changes, please make sure you have [pre-commit](https://pre-commit.com/)
installed and run `pre-commit install` in the module root directory. This installs git hooks
to automatically check coding style and some linters.

## Issues

Take a look at the [issues](https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt/-/issues)
to help shape the future of this project.

# Coding style

Every project has its own style, so here it is: we use clang-format to handle all formatting and provide a
corresponding `.clang-format` style specification file (some of which we explain below). It is convenient
to integrate `clang-format` with your favorite IDE.
This means that everything related to whitespace, brackets, comments, braces, macros, etc. is determined by
clang-format, no flexibilty here (apart from one exception, see below).
_I repeat: no flexibility at all regarding these._
Apart from that we have some simple rules that you are asked to follow and any contributions in violation may be
rejected.
In particular, these are:

* We use `CamelCase` for class and struct names, enum names, typedefs, and the like. We use `stl_standard` for class
methods, members and variables. Does not seem to be consistent, but we find it highly readable. The only exception are typedefs and helper structs which arise in a context similar to the stl, where we use
lowercase `typedef ... type;` and `struct is_vector {`.
* Class members end with a trailing underscore. Static class members may drop the trailing underscore (in particular
dimensions, such as `static const size_t dim = 2;`).
* We append `Type` to in-class typedefs.
* Prefer `using typename BaseClass::FooType;` over `typedef typename BaseType::FooType FooType;`, reduces errors due to
typos.
* Include guards follow a pattern: file location `dune/xt/common/fvector.hh` turns into `DUNE_XT_COMMON_FVECTOR_HH`
* Be careful to handle integers corretly, bad integer conversion is bad! Whenever possible, use `size_t` (even for
`template<class Foo, size_t bar>`) and convert to the correct type using our `numeric_cast` from
[dune/xt/common/numeric_cast.hh](https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt/-/blob/master/dune/xt/common/numeric_cast.hh):
```c++
#include <dune/xt/common/numeric_cast.hh>

const auto unsigned_number = Dune::XT::Common::numeric_cast<size_t>(some_vector_size_of_unknown_int_type);
// work with unsigned_number as size_t
Dune::function_that_requires_an_int_because_ints_are_cool(Dune::XT::Common::numeric_cast<int>(unsigned_number);
```
* Regarding whitespace, the only choice we have is the number of empty lines: one or two.
Within functions and methods, only very rarely place an empty line!
If you feel you need a visual break, placing a comment is often a better idea.
Within a struct or class, place one newline before each method and before the block of members.
Between classes, free functions and the like in a header, always place two newlines to visually separate them.
If you have several template specializations of the same class or function, one newline does suffice.
For example:
```c++
#include "header.hh"

namespace Dune { // one newline between includes and namespaces
namespace XT {
namespace Common {


class Foo // two newlines above this class
{
  typedef Foo ThisType;
public:
  static dim = 1; // this member does not need the trailing underscore

  Foo() // one newline above this ctor
  {}

  void member_func() // one newline above this member function
  {}

private: // one newline above this group of members
  double number_; // all private members with
  bool member_;   // trailing underscore
}; // two newlines below this class


template< int d, bool anything > // two newlines above this free function
void free_function()
{
}

template< int d, bool anything > // but only one above the specialization
void free_function()
{
} // two newlines below this group of specializations


} // namespace Common
} // namespace XT
} // namespace Dune
```
* If you do not have a very good reason against, sort your headers in groups and within a group alphabetically.
Sort the groups by
  * stl
  * other external dependency
  * dune core modules
  * dune-xt
  * local header within the same (sub)folder

For example:
```c++
#include <string>

#include <boost/...>

#include <dune/common/a.hh>
#include <dune/common/b.hh>

#include <dune/xt/common/c.hh>
#include <dune/xt/common/d.hh>

#include "local_a.hh"
#include "local_b.hh"
```


# Project layout

Following the layout of the core modules, all headers are placed below `dune/xt/{common,grid,functions,la}`.
If extending some functionality from dune-common, the names of the headers should coincide, i.e.
`dune/xt/common/fvector.hh` for extensions of `dune/common/fvector.hh`.
In such cases, the first header is expected to include the latter!
Rather use a subdirectory and several headers than overly large headers.

Tests can be found in `dune/xt/test`, where the name of a test source correspond to the full path of the
corresponding header, e.g. `dune/xt/test/common/foo_bar.cc` contains tests for `dune/xt/common/foo/bar.hh`.

