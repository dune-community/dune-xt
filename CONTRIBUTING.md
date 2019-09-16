```
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2017 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
```

# Help

If you are experiencing problems or would like to get help, open up a new
[issue](https://github.com/dune-community/dune-xt/issues/new) and tag it accordingly
or ask for help on the [mailinglist](https://listserv.uni-muenster.de/mailman/listinfo/dune-xt-dev).
Please note that dune-xt-common is not one of the core DUNE modules. Thus, you will not get
any support from the official channels (such as the DUNE mailinglist, the DUNE bugtracker,
etc...).
When submitting bugs, please read these
[general guidelines](https://www.dune-project.org/doc/guides/bug_reporting/) beforehand.


# Contributing

We very much encourage and welcome contributions!
There are several ways to contribute to this project:

## Contribute code and/or documentation

Fork the repository and open up a pull request, explaining your changes and why they should be included.

## Mailinglist and issues

Join the [mailinglist](https://listserv.uni-muenster.de/mailman/listinfo/dune-xt-dev) and take a look at the
[issues](https://github.com/dune-community/dune-xt/issues) to help shape the future of this project.


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
`template<class Foo, size_t bar>`) and convert to the correct type using boost:
```c++
#include <boost/numeric/conversion/cast.hpp>

const auto unsigned_number = boost::numeric_cast<size_t>(some_vector_size_of_unknown_int_type);
// work with unsigned_number as size_t
Dune::function_that_requires_an_int_because_ints_are_cool(boost::numeric_cast<int>(unsigned_number);
```
* Regarding whitespace, the only choice we have is the number of empty lines: one or two.
Within functions and methods, only very rarely place an empty line!
If you feel you need a visual break, placing a comment is often a better idea.
Within a struct or class, place one newline before each method and before the block of members.
Between classes, free functions and the like in a header, always place two newlines to visually seperate them.
If you have several tempalte specializations of the same class or function, one newline does suffice.
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
  double number_; // all provate members with
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
  * dune-xt-common
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

Following the layout of the core modules, all headers are placed below `dune/xt/common`.
If extending some functionality from dune-common, the names of the headers should coincide, i.e.
`dune/xt/common/fvector.hh` for extensions of `dune/common/fvector.hh`.
In such cases, the first header is expected to include the latter!
Rather use a subdirectory and several headers than overly large headers.

Tests can be found in `dune/xt/common/test`, where the name of a test source correspond to the full path of the
corresponding header, e.g. `dune/xt/common/test/foo_bar.cc` contains tests for `dune/xt/common/foo/bar.hh`.


# clang-format

You need to have clang-format in the correct version installed (use `clang-format --version` to check).
To apply clang-format to all files in this project, call `make format` in your build directory.
_Note that clang-format changes the files inplace, so save your work in your editor beforehand!_

A pre-commit hook to test for compliance is automatically installed by ```dunecontrol {vcsetup,all}```.
If you have several versions of clang-format installed you can manually specify the name of the
executble by setting:
```
git config hooks.clangformat my-custom-clang-format-version
```

If you do not have clang-format in the correct version available or forgot to apply it, the code you publish will fail the
corresponding travis check.
Upon failure you are provided with more information on the github pull request page on how to proceed.
In particular you are provided with the required changes which you may download and incorporate into your code.
