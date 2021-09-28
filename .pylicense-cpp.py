# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2019 - 2021)
# ~~~

name = 'This file is part of the dune-xt project:'
url = 'https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt'
copyright_statement = 'Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.'
license = '''Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
          with "runtime exception" (http://www.dune-project.org/license.html)'''
prefix = '//'

include_patterns = ('*.cc', '*.cxx', '*.hh', '*.hxx', '*cmake_config.h.in', '*headercheck.cpp.in', '*config.h.cmake',
                    '*version.hh.in', '*.pbh', '*.tpl')
exclude_patterns = ('*mathexpr.*', '*gtest-all.cc', '*.vcsetup*', '*gtest-all.cxx',
                    '*dune/xt/functions/expression/mathexpr.cc', '*dune/xt/functions/expression/mathexpr.hh')
