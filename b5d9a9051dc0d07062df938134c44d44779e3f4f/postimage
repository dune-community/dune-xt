# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2018)
# ~~~

include(DuneXTTesting)

dune_pybindxi_add_module(dunextusercode EXCLUDE_FROM_ALL dunextusercode.cc)
dune_pybindxi_add_module(dunegdtusercode EXCLUDE_FROM_ALL ${header} dunegdtusercode.cc)

add_subdirectory(dune)
add_subdirectory(scripts)
add_subdirectory(wrapper)

dxt_add_python_tests()
