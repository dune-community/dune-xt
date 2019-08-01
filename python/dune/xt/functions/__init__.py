# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2018)
#   Tim Keil        (2018)
# ~~~

from importlib import import_module

import dune.xt.common
import dune.xt.la
import dune.xt.grid

_modules = (
    '_function_interface_1d',
    '_function_interface_2d',
    '_function_interface_3d',
    '_gridfunction_interface_1d',
    '_gridfunction_interface_2d',
    '_gridfunction_interface_3d',
    '_checkerboard',
    '_constant',
    '_expression',
    '_indicator',
    '_spe10',
    )

# see https://stackoverflow.com/questions/43059267/how-to-do-from-module-import-using-importlib
for mod_name in _modules:
    try:
        mod = import_module('.{}'.format(mod_name), 'dune.xt.functions')
        if "__all__" in mod.__dict__:
            names = mod.__dict__["__all__"]
        else:
            # otherwise we import all names that don't begin with _
            names = [x for x in mod.__dict__ if not x.startswith("_")]
        globals().update({k: getattr(mod, k) for k in names})
    except ImportError as e:
        import os
        import logging
        if os.environ.get('DXT_PYTHON_DEBUG', False):
            raise e
        logging.error('dune-xt-functions: could not import {} module'.format(mod_name))

