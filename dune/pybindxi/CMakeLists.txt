# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2020)
#   Tobias Leibner  (2020 - 2021)
# ~~~

dune_library_add_sources(dunepybindxi SOURCES interpreter.cc)

install(
  FILES attr.h
        buffer_info.h
        cast.h
        chrono.h
        common.h
        complex.h
        eigen.h
        embed.h
        eval.h
        functional.h
        gil.h
        interpreter.hh
        iostream.h
        numpy.h
        operators.h
        options.h
        pybind11.h
        pytypes.h
        stl.h
        stl_bind.h
        stl/filesystem.h
        detail/class.h
        detail/common.h
        detail/descr.h
        detail/init.h
        detail/internals.h
        detail/type_caster_base.h
        detail/typeid.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pybindxi)
