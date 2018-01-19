# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)

from importlib import import_module

from dune.xt.common import DEBUG # inits MPI via mpi4py

_init_logger_methods = list()
_test_logger_methods = list()
_init_mpi_methods = list()
_other_modules = ('xt.common', 'xt.grid', 'xt.functions', 'gdt')

from ._la import __dict__ as module
to_import = [name for name in module if not name.startswith('_')]
globals().update({name: module[name] for name in to_import})
_init_logger_methods.append(module['_init_logger'])
_test_logger_methods.append(module['_test_logger'])
_init_mpi_methods.append(module['_init_mpi'])
del to_import
del module


def init_logger(max_info_level=999,
                max_debug_level=999,
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

def test_logger(info=True, debug=True, warning=True):
    test_logger_methods = _test_logger_methods.copy()
    for module_name in _other_modules:
        try:
            mm = import_module('dune.{}'.format(module_name))
            for test_logger_method in mm._test_logger_methods:
                test_logger_methods.append(test_logger_method)
        except ModuleNotFoundError:
            pass
    for test_logger_method in test_logger_methods:
        test_logger_method(info, debug, warning)

def init_mpi(args=list()):
    if DEBUG:
        init_mpi_methods = [_init_mpi_methods[0],]
    else:
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


HAVE_DUNE_ISTL = 'IstlDenseVectorDouble' in globals()
HAVE_EIGEN = 'EigenDenseVectorDouble' in globals()

