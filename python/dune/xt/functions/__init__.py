# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

import dune.xt

try:
    from dune.xt._interfaces import *
    from dune.xt._functions import *
except ImportError as e:
    import os
    import logging
    if os.environ.get('DXT_PYTHON_DEBUG', False):
        raise e
    logging.error('dune-xt-functions bindings not available')
