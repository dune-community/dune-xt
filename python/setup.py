#!/usr/bin/env python
#
# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Rene Fritze (2018)
#
# This file is part of the dune-xt-common project:
# ~~~

import sys
from setuptools import setup, find_packages

setup(name='dune.xt.la',
      version='2.4',
      namespace_packages=['dune'],
      description='Python for Dune-Xt',
      author='The dune-xt devs',
      author_email='dune-xt-dev@listserv.uni-muenster.de',
      url='https://github.com/dune-community/dune-xt-la',
      install_requires=['numpy', 'scipy'],
      packages = find_packages(),
      zip_safe = 0,
      package_data = {'': ['*.so']},)
