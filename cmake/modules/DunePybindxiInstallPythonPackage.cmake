# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2018)
#   René Fritze     (2018, 2020)
#   Tobias Leibner  (2018, 2020)
# ~~~

macro(INCLUDE_DEPENDENT_BINARY_PYTHON_DIRS) # disable most warnings from dependent modules
  foreach(mod ${ALL_DEPENDENCIES} ${PROJECT_NAME})
    dune_module_path(MODULE ${mod} RESULT ${mod}_binary_dir BUILD_DIR)
    set(tdir ${${mod}_binary_dir})
    if(IS_DIRECTORY ${tdir})
      dune_register_package_flags(INCLUDE_DIRS ${tdir})
    endif()
  endforeach(mod DEPENDENCIES)
endmacro(INCLUDE_DEPENDENT_BINARY_PYTHON_DIRS)

# Copy of dune_python_install_package from dune-common. Changes:
#
# * Package is always installed into the dune-env, even if a setup.py.in is found instead of a setup.py.
# * If a setup.py.in is found, the whole directory is symlinked to the binary dir.
#
# cmake-lint: disable=R0915
function(dune_pybindxi_install_python_package)
  # Parse arguments
  set(option)
  set(single PATH)
  set(multi ADDITIONAL_PIP_PARAMS)
  include(CMakeParseArguments)
  cmake_parse_arguments(PYINST "${option}" "${single}" "${multi}" ${ARGN})
  if(PYINST_UNPARSED_ARGUMENTS)
    message(WARNING "Unparsed arguments in dune_python_install_package: This often indicates typos!")
  endif()

  # Check for the presence of the pip package
  dune_python_find_package(PACKAGE pip INTERPRETER ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE})
  if(NOT DUNE_PYTHON_pip_FOUND)
    message(FATAL_ERROR "dune_python_install_package: Requested installations, but pip was not found!")
  endif()

  set(pyinst_fullpath ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH})
  if(EXISTS ${pyinst_fullpath}/setup.py.in)
    # symlink files to binary dir
    file(GLOB_RECURSE files "${pyinst_fullpath}/*")
    foreach(fn ${files})
      file(RELATIVE_PATH rel_fn ${CMAKE_CURRENT_SOURCE_DIR} ${fn})
      get_filename_component(directory ${rel_fn} DIRECTORY)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${directory})
      execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${rel_fn}"
                              "${CMAKE_CURRENT_BINARY_DIR}/${rel_fn}")
    endforeach()
    configure_file(${PYINST_PATH}/setup.py.in ${PYINST_PATH}/setup.py)
    set(pyinst_fullpath ${CMAKE_CURRENT_BINARY_DIR}/${PYINST_PATH})
    set(pyinst_purepython TRUE)
  elseif(EXISTS ${pyinst_fullpath}/setup.py)
    set(pyinst_purepython TRUE)
  else()
    message(
      FATAL_ERROR "dune_python_install_package: Requested installations, but neither setup.py nor setup.py.in found!")
  endif()

  # Find out whether we should install in editable mode
  set(install_editable 1)

  # Construct the wheel house installation option string
  set(wheel_option "")
  if(IS_DIRECTORY ${DUNE_PYTHON_WHEELHOUSE})
    set(wheel_option "--find-links=${DUNE_PYTHON_WHEELHOUSE}")
    #
    # The following line is a bummer! We cannot have editable packages once we start using global installations! This is
    # related to the nightmare that is https://github.com/pypa/pip/issues/3
    #
    set(install_editable FALSE)
  endif()

  # Construct the editable option string
  set(edit_option "")
  if(install_editable)
    set(edit_option "-e")
  endif()

  # Construct the installation location option string
  set(install_option "")
  if("${DUNE_PYTHON_INSTALL_LOCATION}" STREQUAL "user")
    set(install_option "--user")
  endif()

  set(install_cmdline
      -m
      pip
      install
      "${install_option}"
      "${wheel_option}"
      "${edit_option}"
      ${PYINST_ADDITIONAL_PIP_PARAMS}
      "${pyinst_fullpath}")

  #
  # If requested, install into the configure-time Dune virtualenv
  #

  if(pyinst_purepython AND DUNE_PYTHON_VIRTUALENV_SETUP)
    message("-- Installing python package at ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH} into the virtualenv...")
    dune_execute_process(COMMAND "${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}" "${install_cmdline}" ERROR_MESSAGE
                         "dune_python_install_package: Error installing into virtualenv!")
  endif()

  #
  # Now define rules for `make install_python`.
  #

  # Leave this function if no installation rules are required
  dune_module_path(MODULE dune-common RESULT scriptdir SCRIPT_DIR)

  # Determine a target name for installing this package
  string(REPLACE "/" "_" targetname "install_python_${CMAKE_CURRENT_SOURCE_DIR}_${PYINST_PATH}")

  # Add a custom target that globally installs this package if requested
  add_custom_target(
    ${targetname}
    COMMAND ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} ${install_cmdline}
    COMMENT "Installing the python package at ${pyinst_fullpath}")

  add_dependencies(install_python ${targetname})

  # Define rules for `make install` that install a wheel into a central wheelhouse
  #
  # NB: This is necessary, to allow mixing installed and non-installed modules with python packages. The wheelhouse will
  # allow to install any missing python packages into a virtual environment.
  #

  # Construct the wheel installation commandline
  set(wheel_command ${PYTHON_EXECUTABLE} -m pip wheel -w ${DUNE_PYTHON_WHEELHOUSE} ${pyinst_fullpath})

  # Add the installation rule
  install(
    CODE "message(\"Installing wheel for python package at ${pyinst_fullpath}...\")
                dune_execute_process(COMMAND ${wheel_command}
                                     ERROR_MESSAGE \"Error installing wheel for python package at ${pyinst_fullpath}\"
                                     )")
endfunction()
