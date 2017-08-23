#!/usr/bin/env python
#
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2017 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2016 - 2017)

import sys
from setuptools import setup

setup(name='dune.xt.functions',
      version='2.4',
      namespace_packages=['dune'],
      description='Python for Dune-Xt',
      author='The dune-xt devs',
      author_email='dune-xt-dev@listserv.uni-muenster.de',
      url='https://github.com/dune-community/dune-xt-functions',
      packages=['dune.xt'],
      install_requires=[],
      scripts=[])
