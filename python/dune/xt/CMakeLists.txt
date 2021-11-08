# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019 - 2020)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018)
#   Tobias Leibner  (2019 - 2020)
# ~~~

add_subdirectory(common)
add_subdirectory(functions)
add_subdirectory(grid)
add_subdirectory(la)
add_subdirectory(test)

configure_file(_version.py.in _version.py)
