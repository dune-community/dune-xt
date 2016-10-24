# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

add_subdirectory(la)

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(la la/bindings.cc)
  target_link_dune_default_libraries(la)
endif()

