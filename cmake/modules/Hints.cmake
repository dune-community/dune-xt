# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018 - 2020)
# ~~~

macro(append_to_each INPUTLIST POSTFIX OUTPUTLIST)
  foreach(ENTRY ${INPUTLIST})
    list(APPEND ${OUTPUTLIST} ${ENTRY}${POSTFIX})
  endforeach(ENTRY ${INPUTLIST})
endmacro()

set(root_hints
    "${CMAKE_SOURCE_DIR}/../local/"
    "${CMAKE_SOURCE_DIR}/../environments/debian-minimal/local/"
    "${CMAKE_SOURCE_DIR}/../environments/debian-full/local/"
    "$ENV{HOME}/"
    "$ENV{HOME}/Software/")

set(bin_hints "")
append_to_each("${root_hints}" "bin/" bin_hints)

set(lib_hint "")
append_to_each("${root_hints}" "lib/" lib_hints)

set(include_hints "")
append_to_each("${root_hints}" "include/" include_hints)
