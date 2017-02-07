# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

set(lib_dune_xt_la_sources
    container/pattern.cc)

dune_library_add_sources(dunextla SOURCES ${lib_dune_xt_la_sources})

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(_la bindings.cc)
  target_link_dune_default_libraries(_la)
  if(DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS _la)
  else()
    add_custom_target(bindings DEPENDS _la)
  endif()
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)


