# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2016, 2018)
# ~~~

# Note: This import makes sure that metis is imported before python bindings from dune-xt. Importing dune-xt python
# bindings at any point before metis will cause a segfault due to misscommunication between libmetis5 and
# libscotchmetis.

try:
	import metis
except ImportError:
	pass


__import__('pkg_resources').declare_namespace(__name__)


from importlib import import_module
import os
import logging


def guarded_import(globs, base_name, mod_name):
    # see https://stackoverflow.com/questions/43059267/how-to-do-from-module-import-using-importlib
    try:
        mod = import_module('.{}'.format(mod_name), base_name)
        if "__all__" in mod.__dict__:
            names = mod.__dict__["__all__"]
        else:
            # otherwise we import all names that don't begin with _
            names = [x for x in mod.__dict__ if not x.startswith("_")]
        for nm in names:
            if nm in globs:
                logging.error('{}: overwriting existing name \'{}\' when importing from \'{}\' (continuing
                        anyway)!'.format(base_name, nm, mod_name))
        globs.update({k: getattr(mod, k) for k in names})
    except ImportError as e:
        logging.error('{}: could not import module \'{}\' (continuing anyway)!'.format(base_name, mod_name))
        if os.environ.get('DXT_PYTHON_DEBUG', False):
            raise e

