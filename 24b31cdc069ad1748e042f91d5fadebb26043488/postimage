# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2018)
#   René Fritze     (2018)
#   Tobias Leibner  (2018)
# ~~~

add_subdirectory(xt)
add_subdirectory(pybindxi)
install(
  DIRECTORY xt
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/
  FILES_MATCHING
  PATTERN "*.hh")
add_subdirectory(gdt)
