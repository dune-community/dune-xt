# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2018)
# ~~~

name = 'This file is part of the dune-xt-grid project:'
url = 'https://github.com/dune-community/dune-xt-grid'
copyright_statement = 'Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.'
license = '''Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
          with "runtime exception" (http://www.dune-project.org/license.html)'''
prefix = '#'
lead_in = '# ~~~'
lead_out = '# ~~~'

include_patterns = ('*.txt', '*.cmake', '*.py', '*.sh', '*.bash', '*.dgf', '*.msh', '*.gdb', '*.cfg',
                    '*.travis.*', '*.gitignore', '*.mailmap', '*.gitattributes', '*gitignore-*', '*stamp-vc',
                    '*dune.module', '*Doxylocal', '*.clang-format', '*COPYING-CMAKE-SCRIPTS', '*README',
                    '*LICENSE', '*mainpage', '*switch-build_dir', '*dune-xt-grid.pc.in', '*CMakeLists.txt')
exclude_patterns = ('*config.h.cmake', '*.vcsetup*', '*.dgf', '*builder_definitions.cmake')

