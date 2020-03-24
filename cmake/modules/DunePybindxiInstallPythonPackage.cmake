# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
#
# ~~~

macro(include_dependent_binary_python_dirs) # disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES} ${PROJECT_NAME})
    dune_module_path(MODULE ${_mod} RESULT ${_mod}_binary_dir BUILD_DIR)
    set(tdir ${${_mod}_binary_dir})
    if(IS_DIRECTORY ${tdir})
      dune_register_package_flags(INCLUDE_DIRS ${tdir})
    endif()
  endforeach(_mod DEPENDENCIES)
endmacro(include_dependent_binary_python_dirs)

# copy of dune_python_install_package from dune-common
# changes: - package is always installed into the dune-env, even if a setup.py.in is found instead of a setup.py
#          - if a setup.py.in is found, the whole directory is symlinked to the binary dir
function(dune_pybindxi_install_python_package)
  # Parse Arguments
  set(OPTION)
  set(SINGLE PATH)
  set(MULTI ADDITIONAL_PIP_PARAMS)
  include(CMakeParseArguments)
  cmake_parse_arguments(PYINST "${OPTION}" "${SINGLE}" "${MULTI}" ${ARGN})
  if(PYINST_UNPARSED_ARGUMENTS)
    message(WARNING "Unparsed arguments in dune_python_install_package: This often indicates typos!")
  endif()

  # Check for the presence of the pip package
  if(NOT DUNE_PYTHON_pip_FOUND)
    message(FATAL_ERROR "dune_python_install_package: Requested installations, but pip was not found!")
  endif()

  set(PYINST_FULLPATH ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH})
  if(EXISTS ${PYINST_FULLPATH}/setup.py.in)
    # symlink files to binary dir
    file(GLOB_RECURSE files "${PYINST_FULLPATH}/*")
    foreach(fn ${files})
      file(RELATIVE_PATH rel_fn ${CMAKE_CURRENT_SOURCE_DIR} ${fn})
      get_filename_component(directory ${rel_fn} DIRECTORY)
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${directory})
      execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${CMAKE_CURRENT_SOURCE_DIR}/${rel_fn}"
	      "${CMAKE_CURRENT_BINARY_DIR}/${rel_fn}")
    endforeach()
    configure_file(${PYINST_PATH}/setup.py.in ${PYINST_PATH}/setup.py)
    set(PYINST_FULLPATH ${CMAKE_CURRENT_BINARY_DIR}/${PYINST_PATH})
    set(PYINST_PUREPYTHON TRUE)
  elseif(EXISTS ${PYINST_FULLPATH}/setup.py)
    set(PYINST_PUREPYTHON TRUE)
  else()
    message(FATAL_ERROR "dune_python_install_package: Requested installations, but neither setup.py nor setup.py.in found!")
  endif()

  # Find out whether we should install in editable mode
  set(INSTALL_EDITABLE 1)

  # Construct the wheel house installation option string
  set(WHEEL_OPTION "")
  if(IS_DIRECTORY ${DUNE_PYTHON_WHEELHOUSE})
    set(WHEEL_OPTION "--find-links=${DUNE_PYTHON_WHEELHOUSE}")
    #
    # The following line is a bummer!
    # We cannot have editable packages once we start using global installations!
    # This is related to the nightmare that is https://github.com/pypa/pip/issues/3
    #
    set(INSTALL_EDITABLE FALSE)
  endif()

  # Construct the editable option string
  set(EDIT_OPTION "")
  if(INSTALL_EDITABLE)
    set(EDIT_OPTION "-e")
  endif()

  # Construct the installation location option string
  set(INSTALL_OPTION "")
  if("${DUNE_PYTHON_INSTALL_LOCATION}" STREQUAL "user")
    set(INSTALL_OPTION "--user")
  endif()

  set(INSTALL_CMDLINE -m pip install
                      "${INSTALL_OPTION}" "${WHEEL_OPTION}" "${EDIT_OPTION}" ${PYINST_ADDITIONAL_PIP_PARAMS}
                      "${PYINST_FULLPATH}")


  #
  # If requested, install into the configure-time Dune virtualenv
  #

  if(PYINST_PUREPYTHON AND DUNE_PYTHON_VIRTUALENV_SETUP)
    message("-- Installing python package at ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH} into the virtualenv...")
    dune_execute_process(COMMAND "${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}" "${INSTALL_CMDLINE}"
                         ERROR_MESSAGE "dune_python_install_package: Error installing into virtualenv!")
  endif()

  #
  # Now define rules for `make install_python`.
  #

  # Leave this function if no installation rules are required
  dune_module_path(MODULE dune-common
                   RESULT scriptdir
                   SCRIPT_DIR)

  # Determine a target name for installing this package
  string(REPLACE "/" "_" targetname "install_python_${CMAKE_CURRENT_SOURCE_DIR}_${PYINST_PATH}")

  # Add a custom target that globally installs this package if requested
  add_custom_target(${targetname}
	            COMMAND ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} ${INSTALL_CMDLINE}
                    COMMENT "Installing the python package at ${PYINST_FULLPATH}"
                    )

  add_dependencies(install_python ${targetname})

  # Define rules for `make install` that install a wheel into a central wheelhouse
  #
  # NB: This is necessary, to allow mixing installed and non-installed modules
  #     with python packages. The wheelhouse will allow to install any missing
  #     python packages into a virtual environment.
  #

  # Construct the wheel installation commandline
  set(WHEEL_COMMAND ${PYTHON_EXECUTABLE} -m pip wheel -w ${DUNE_PYTHON_WHEELHOUSE} ${PYINST_FULLPATH})

  # Add the installation rule
  install(CODE "message(\"Installing wheel for python package at ${PYINST_FULLPATH}...\")
                dune_execute_process(COMMAND ${WHEEL_COMMAND}
                                     ERROR_MESSAGE \"Error installing wheel for python package at ${PYINST_FULLPATH}\"
                                     )"
          )
endfunction()
