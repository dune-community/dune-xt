# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

install(FILES
        attr.h
        buffer_info.h
        cast.h
        chrono.h
        common.h
        complex.h
        eigen.h
        embed.h
        eval.h
        functional.h
        iostream.h
        numpy.h
        operators.h
        options.h
        pybind11.h
        pytypes.h
        stl.h
        stl_bind.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/pybindxi)

