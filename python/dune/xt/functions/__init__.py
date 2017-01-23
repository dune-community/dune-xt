# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# The copyright lies with the authors of this file (see below).
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)

from importlib import import_module


def init_logger(max_info_level=-1,
                max_debug_level=-1,
                enable_warnings=True,
                enable_colors=True,
                info_color='blue',
                debug_color='darkgray',
                warning_color='red'):
    from ._functions import init_logger as _init_logger
    initializers = [_init_logger]
    for module_name in ('xt.common', 'xt.grid', 'xt.la', 'gdt'):
        try:
            mm = import_module('dune.{}'.format(module_name))
            initializers.append(mm.init_logger)
        except ModuleNotFoundError:
            pass
    for initializer in initializers:
        initializer(max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color,
                        warning_color)


def init_mpi(args=list()):
    try:
        from dune.gdt import init_mpi as _init_mpi
    except ModuleNotFoundError:
        from dune.xt.common import init_mpi as _init_mpi
    _init_mpi(args)


from ._functions import *
