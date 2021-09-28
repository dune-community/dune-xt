# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Tobias Leibner  (2021)
# ~~~

set(DXT_DATA_BASEDIR "${CMAKE_BINARY_DIR}/data")
if (Spe10Data_FOUND)
  return()
endif()
set(Spe10Data_FOUND 0)
if (NOT EXISTS "${DXT_DATA_BASEDIR}")
  file(MAKE_DIRECTORY "${DXT_DATA_BASEDIR}")
endif(NOT EXISTS "${DXT_DATA_BASEDIR}")
if (NOT EXISTS "${DXT_DATA_BASEDIR}/perm_case1.dat")
  # download file 1
  file(DOWNLOAD https://github.com/dune-community/dune-xt-data/raw/master/dune/xt/data/perm_case1.zip ${DXT_DATA_BASEDIR}/perm_case1.zip STATUS DOWNLOAD_STATUS)
  list(GET DOWNLOAD_STATUS 0 SPE10_FILE1_STATUS)
  if (NOT SPE10_FILE1_STATUS EQUAL 0)
    message(WARNING "Failed downloading perm_case1.zip! Error: ${DOWNLOAD_STATUS}.")
    set(Spe10Data_FOUND 0)
    return()
  endif(NOT SPE10_FILE1_STATUS EQUAL 0)
  # extract file 1
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${DXT_DATA_BASEDIR}/perm_case1.zip"
    WORKING_DIRECTORY "${DXT_DATA_BASEDIR}"
    RESULT_VARIABLE SPE10_FILE1_STATUS
    OUTPUT_VARIABLE SPE10_FILE1_OUTPUT
    ERROR_VARIABLE SPE10_FILE1_OUTPUT)
  if (NOT SPE10_FILE1_STATUS EQUAL 0)
    message(WARNING "Failed extracting perm_case1.zip! Output: ${SPE10_FILE1_OUTPUT}.")
    set(Spe10Data_FOUND 0)
    return()
  endif(NOT SPE10_FILE1_STATUS EQUAL 0)
endif(NOT EXISTS ${DXT_DATA_BASEDIR}/perm_case1.dat)

if (NOT EXISTS "${DXT_DATA_BASEDIR}/spe_perm.dat")
  # download file 2
  file(DOWNLOAD https://github.com/dune-community/dune-xt-data/raw/master/dune/xt/data/por_perm_case2a.zip ${DXT_DATA_BASEDIR}/por_perm_case2a.zip STATUS DOWNLOAD_STATUS)
  list(GET DOWNLOAD_STATUS 0 SPE10_FILE2_STATUS)
  if (NOT SPE10_FILE2_STATUS EQUAL 0)
    message(WARNING "Failed downloading por_perm_case2a.zip! Error: ${DOWNLOAD_STATUS}.")
    set(Spe10Data_FOUND 0)
    return()
  endif(NOT SPE10_FILE2_STATUS EQUAL 0)
  # extract file 2
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${DXT_DATA_BASEDIR}/por_perm_case2a.zip"
    WORKING_DIRECTORY "${DXT_DATA_BASEDIR}"
    RESULT_VARIABLE SPE10_FILE2_STATUS
    OUTPUT_VARIABLE SPE10_FILE2_OUTPUT
    ERROR_VARIABLE SPE10_FILE2_OUTPUT)
  if (NOT SPE10_FILE2_STATUS EQUAL 0)
    message(WARNING "Failed extracting por_perm_case2a.zip! Output: ${SPE10_FILE2_OUTPUT}.")
    set(Spe10Data_FOUND 0)
    return()
  endif(NOT SPE10_FILE2_STATUS EQUAL 0)
endif (NOT EXISTS "${DXT_DATA_BASEDIR}/spe_perm.dat")
set(Spe10Data_FOUND 1)
set(HAVE_SPE10_DATA 1)
set(SPE10_MODEL1_FILENAME ${DXT_DATA_BASEDIR}/perm_case1.dat)
set(SPE10_MODEL2_FILENAME ${DXT_DATA_BASEDIR}/spe_perm.dat)

