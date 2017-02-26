# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

set(lib_dune_xt_functions_sources
    expression/mathexpr.cc)

dune_library_add_sources(dunextfunctions SOURCES ${lib_dune_xt_functions_sources})

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(_functions bindings.cc)
  target_link_dune_default_libraries(_functions)
  add_dune_alberta_flags(GRIDDIM 2 _functions)
  if (DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS _functions)
  else()
    add_custom_target(bindings DEPENDS _functions)
  endif()
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)

