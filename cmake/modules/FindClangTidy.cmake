# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2016, 2018 - 2019)
#
# Use this module by invoking find_package with the form::
#
#   find_package(ClangTidy
#     [version] # Minimum version e.g. 3.7
#     )
#
#   optionally pass the minimal version you require like so find_package(ClangTidy 3.7)
#   this module set ClangTidy_EXECUTABLE, ClangTidy_VERSION
#   and ClangTidy_FOUND accordingly
# ~~~

find_program(ClangTidy_EXECUTABLE
             NAMES clang-tidy
                   clang-tidy-3.6
                   clang-tidy-3.7
                   clang-tidy-3.8
                   clang-tidy-3.9
                   clang-tidy-4.0)
if(EXISTS ${ClangTidy_EXECUTABLE})
  execute_process(COMMAND ${ClangTidy_EXECUTABLE} -version OUTPUT_VARIABLE clang_out)
  string(REGEX
         REPLACE ".*LLVM version ([0-9]+\\.[0-9]+).*"
                 "\\1"
                 ClangTidy_VERSION
                 ${clang_out})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ClangTidy REQUIRED_VARS ClangTidy_EXECUTABLE VERSION_VAR ClangTidy_VERSION)
