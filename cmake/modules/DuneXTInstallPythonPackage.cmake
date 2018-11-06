# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
#
# ~~~

# copy from dune-common
function(dune_xt_execute_process)
  include(CMakeParseArguments)
  cmake_parse_arguments(EXECUTE "" "ERROR_MESSAGE;RESULT_VARIABLE;OUTPUT_VARIABLE;ERROR_VARIABLE" "" ${ARGN})

  # Decide whether stdout and stderr have to be split
  if(EXECUTE_OUTPUT_VARIABLE AND EXECUTE_ERROR_VARIABLE)
    set(SPLIT_ERROR TRUE)
    set(ERRLOGVAR errlog)
  else()
    set(SPLIT_ERROR FALSE)
    set(ERRLOGVAR log)
  endif()

  # Call the original cmake function
  execute_process(${EXECUTE_UNPARSED_ARGUMENTS}
                  RESULT_VARIABLE retcode
                  OUTPUT_VARIABLE log
                  ERROR_VARIABLE ${ERRLOGVAR}
                  )

  # Issue an error if requested!
  if(EXECUTE_ERROR_MESSAGE)
    if(NOT "${retcode}" STREQUAL "0")
      cmake_parse_arguments(ERR "" "" "COMMAND" ${EXECUTE_UNPARSED_ARGUMENTS})
      if(SPLIT_ERROR)
        set(log "stdout:\n${log}\n\nstderr:\b${errlog}")
      endif()
      message(FATAL_ERROR "${EXECUTE_ERROR_MESSAGE}\nRun command:${ERR_COMMAND}\nReturn code: ${retcode}\nDetailed log:\n${log}")
    endif()
  endif()

  # Propagate logs back to the calling scope
  if(EXECUTE_RESULT_VARIABLE)
    set(${EXECUTE_RESULT_VARIABLE} ${retcode} PARENT_SCOPE)
  endif()
  if(EXECUTE_OUTPUT_VARIABLE)
    set(${EXECUTE_OUTPUT_VARIABLE} ${log} PARENT_SCOPE)
  endif()
  if(EXECUTE_ERROR_VARIABLE)
    set(${EXECUTE_ERROR_VARIABLE} ${${ERROR_VARIABLE}} PARENT_SCOPE)
  endif()
endfunction()

macro(include_dependent_binary_python_dirs) # disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES} ${PROJECT_NAME})
    dune_module_path(MODULE ${_mod} RESULT ${_mod}_binary_dir BUILD_DIR)
    set(tdir ${${_mod}_binary_dir})
    if(IS_DIRECTORY ${tdir})
      dune_register_package_flags(INCLUDE_DIRS ${tdir})
    endif()
  endforeach(_mod DEPENDENCIES)
endmacro(include_dependent_binary_python_dirs)

function(dune_xt_install_python_package) # Parse Arguments
  set(OPTION NO_PIP NO_EDIT)
  set(SINGLE PATH)
  set(MULTI ADDITIONAL_PIP_PARAMS)
  include(CMakeParseArguments)
  cmake_parse_arguments(PYINST "${OPTION}" "${SINGLE}" "${MULTI}" ${ARGN})
  if(PYINST_UNPARSED_ARGUMENTS)
    message(WARNING "Unparsed arguments in dune_install_python_package: This often indicates typos!")
  endif()

  # Determine a target name for installing this package
  string(REPLACE "/"
                 "_"
                 name_suffix
                 ${PYINST_PATH})
  set(targetname "install_python_${DUNE_MOD_NAME}_${name_suffix}")
  if(TARGET ${targetname})
    return()
  endif()

  file(GLOB_RECURSE files "${PROJECT_SOURCE_DIR}/${PYINST_PATH}/*")
  if(PROJECT_SOURCE_DIR EQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "trying to configure with src==bindir")
  endif()
  foreach(fn ${files})
    file(RELATIVE_PATH rel_fn ${PROJECT_SOURCE_DIR} ${fn})
    get_filename_component(directory ${rel_fn} DIRECTORY)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${directory})
    execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${PROJECT_SOURCE_DIR}/${rel_fn}"
                            "${PROJECT_BINARY_DIR}/${rel_fn}")
  endforeach()
  set(PYINST_PATH "${CMAKE_CURRENT_BINARY_DIR}/${PYINST_PATH}")
  dune_register_package_flags(INCLUDE_DIRS ${PYINST_PATH})

  #
  # Install the given python package into dune-python's virtualenv
  #

  # Construct the installation command strings from the given options
  set(WHEEL_ARG "")
  if(IS_DIRECTORY ${DUNE_PYTHON_WHEELHOUSE})
    set(WHEEL_ARG "--find-links=${DUNE_PYTHON_WHEELHOUSE}") # The following line is a bummer!  We cannot have editable
                                                            # packages once we start using global installations!  This
                                                            # is related to the nightmare that is
                                                            # https://github.com/pypa/pip/issues/3
    set(PYINST_NO_EDIT TRUE)
  endif()
  if(PYINST_NO_PIP)
    if(PYINST_NO_EDIT)
      set(INST_COMMAND install)
    else()
      set(INST_COMMAND develop)
    endif()
    set(VENV_INSTALL_COMMAND setup.py ${INST_COMMAND} ${WHEEL_ARG})
  else()
    set(EDIT_OPTION)
    if(NOT PYINST_NO_EDIT)
      set(EDIT_OPTION -e)
    endif()
    set(VENV_INSTALL_COMMAND -m pip install ${PYINST_ADDITIONAL_PIP_PARAMS} ${EDIT_OPTION} ${WHEEL_ARG} .)
  endif()

  # install the package into the virtual env
  dune_xt_execute_process(COMMAND
	  ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}
                          ${VENV_INSTALL_COMMAND}
                          WORKING_DIRECTORY
                          ${PYINST_PATH}
                          ERROR_MESSAGE
                          "Fatal error when installing the package at ${PYINST_PATH} into the env.")

  # Construct the interpreter options for global installation
  if(PYINST_NO_PIP)
    set(SYSTEM_INSTALL_CMDLINE ${PYTHON_EXECUTABLE} setup.py install ${WHEEL_ARG})
    if(DUNE_PYTHON_INSTALL_USER)
      message("Error message: Incompatible options - NO_PIP and DUNE_PYTHON_INSTALL_USER")
    endif()
  else()
    set(USER_STRING "")
    if(DUNE_PYTHON_INSTALL_USER)
      set(USER_STRING --user)
    endif()
    set(SYSTEM_INSTALL_CMDLINE
        ${PYTHON_EXECUTABLE}
        -m
        pip
        install
        ${USER_STRING}
        ${PYINST_ADDITIONAL_PIP_PARAMS}
        ${WHEEL_ARG}
        ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH})
  endif()

  #
  # Now define rules for `make install_python`.
  #

  dune_module_path(MODULE dune-python RESULT DUNE_PYTHON_MODULE_DIR CMAKE_MODULES)

  # Add a custom target that globally installs this package if requested
  add_custom_target(${targetname}
                    COMMAND ${CMAKE_COMMAND} -DCMAKE_MODULE_PATH=${DUNE_PYTHON_MODULE_DIR}
                            -DCMDLINE="${SYSTEM_INSTALL_CMDLINE}" -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} -P
                            ${DUNE_PYTHON_MODULE_DIR}/install_python_package.cmake
                    COMMENT "Installing the python package at ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH}")

  add_dependencies(install_python ${targetname})

  #
  # Define rules for `make install` that install a wheel into a central wheelhouse
  #
  # NB: This is necessary, to allow mixing installed and non-installed modules with python packages. The wheelhouse
  # will allow to install any missing python packages into the virtualenv.
  #

  # Construct the wheel installation commandline
  if(PYINST_NO_PIP)
    set(WHEEL_COMMAND setup.py bdist_wheel)
  else()
    set(WHEEL_COMMAND -m pip wheel)
  endif()

  set(WHEEL_COMMAND
      ${PYTHON_EXECUTABLE}
      ${WHEEL_COMMAND}
      -w
      ${DUNE_PYTHON_WHEELHOUSE}
      ${WHEEL_ARG}
      ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH})

  install(
    CODE
    "message(\"Installing wheel for python package at ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH}...\")
                dune_execute_process(COMMAND ${WHEEL_COMMAND}
                                     ERROR_MESSAGE \"Error installing wheel for python package at ${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH}\"
                                     )")

  #
  # Set some paths needed for Sphinx documentation.
  #

  # Use a custom section to export python path to downstream modules
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
  set(DUNE_PYTHON_SOURCE_PATHS \"${DUNE_PYTHON_SOURCE_PATHS}:${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH}\")
  " PARENT_SCOPE)

  # and add python path for this module
  set(DUNE_PYTHON_SOURCE_PATHS "${DUNE_PYTHON_SOURCE_PATHS}:${CMAKE_CURRENT_SOURCE_DIR}/${PYINST_PATH}" PARENT_SCOPE)

endfunction()
