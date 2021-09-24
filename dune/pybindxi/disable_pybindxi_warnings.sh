# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Tobias Leibner (2020 - 2021)
#
# specifies all pybind11 headers as system headers to disable warnings
# should also work for clang which supports some GCC macros for compatibility
# including as system folder does not work as we would have to include the whole dune-xt folder using -isystem
# ~~~

sed -i '1 i\#pragma GCC system_header' *.h
sed -i '1 i\#pragma GCC system_header' detail/*.h
sed -i '1 i\#pragma GCC system_header' stl/*.h
