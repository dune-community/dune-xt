# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

set(lib_dune_xt_functions_sources
    expression/mathexpr.cc)

dune_library_add_sources(dunextfunctions SOURCES ${lib_dune_xt_functions_sources})

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(_functions bindings.cc)
  target_link_dune_default_libraries(_functions)
  if (DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS _functions)
  else()
    add_custom_target(bindings DEPENDS _functions)
  endif()
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)

