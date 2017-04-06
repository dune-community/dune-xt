# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)

from importlib import import_module

import dune.xt.common # inits MPI via mpi4py

_init_logger_methods = list()
_init_mpi_methods = list()
_other_modules = ('xt.common', 'xt.la', 'xt.functions', 'gdt')

from ._grid import __dict__ as module
to_import = [name for name in module if not name.startswith('_')]
globals().update({name: module[name] for name in to_import})
_init_logger_methods.append(module['_init_logger'])
_init_mpi_methods.append(module['_init_mpi'])
del to_import
del module


def init_logger(max_info_level=-1,
                max_debug_level=-1,
                enable_warnings=True,
                enable_colors=True,
                info_color='blue',
                debug_color='darkgray',
                warning_color='red'):
    init_logger_methods = _init_logger_methods.copy()
    for module_name in _other_modules:
        try:
            mm = import_module('dune.{}'.format(module_name))
            for init_logger_method in mm._init_logger_methods:
                init_logger_methods.append(init_logger_method)
        except ModuleNotFoundError:
            pass
    for init_logger_method in init_logger_methods:
        init_logger_method(max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color,
                           warning_color)

def init_mpi(args=list()):
    init_mpi_methods = _init_mpi_methods.copy()
    for module_name in _other_modules:
        try:
            mm = import_module('dune.{}'.format(module_name))
            for init_mpi_method in mm._init_mpi_methods:
                init_mpi_methods.append(init_mpi_method)
        except ModuleNotFoundError:
            pass
    for init_mpi_method in init_mpi_methods:
        init_mpi_method(args)

