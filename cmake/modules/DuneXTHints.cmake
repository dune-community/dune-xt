# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018 - 2020)
# ~~~

# Appends postfix to each element in inputlist
macro(APPEND_TO_EACH inputlist postfix outputlist)
  foreach(entry ${inputlist})
    list(APPEND ${outputlist} ${entry}${postfix})
  endforeach(entry ${inputlist})
endmacro()

set(_root_hints "/usr/" "${CMAKE_SOURCE_DIR}/../local/" "${CMAKE_SOURCE_DIR}/../environments/debian-minimal/local/"
                "${CMAKE_SOURCE_DIR}/../environments/debian-full/local/" "$ENV{HOME}/" "$ENV{HOME}/Software/")

set(BIN_HINTS "")
append_to_each("${_root_hints}" "bin/" BIN_HINTS)

set(LIB_HINT "")
append_to_each("${_root_hints}" "lib/" LIB_HINTS)

set(INCLUDE_HINTS "")
append_to_each("${_root_hints}" "include/" INCLUDE_HINTS)
