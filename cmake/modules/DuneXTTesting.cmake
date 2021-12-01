# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2017, 2019 - 2020)
#   René Fritze     (2010 - 2020)
#   Sven Kaulmann   (2013)
#   Tobias Leibner  (2015 - 2020)
# ~~~

macro(DXT_HEADERCHECK_TARGET_NAME arg)
  string(REGEX REPLACE ".*/([^/]*)" "\\1" simple ${arg})
  string(REPLACE ${PROJECT_SOURCE_DIR} "" rel ${arg})
  string(REGEX REPLACE "(.*)/[^/]*" "\\1" relpath ${rel})
  string(REGEX REPLACE "/" "_" targname ${rel})
  set(targname "headercheck_${targname}")
endmacro(DXT_HEADERCHECK_TARGET_NAME)

macro(GET_HEADERCHECK_TARGETS subdir)
  file(GLOB_RECURSE bindir_header "${CMAKE_BINARY_DIR}/*.hh")
  list(APPEND dxt_ignore_header ${bindir_header})

  if(ENABLE_HEADERCHECK)
    file(GLOB_RECURSE headerlist "${CMAKE_SOURCE_DIR}/dune/xt/${subdir}/*.hh"
         "${CMAKE_SOURCE_DIR}/dune/xt/test/${subdir}/*.hh" "${CMAKE_SOURCE_DIR}/python/dune/xt/${subdir}/*.hh")
    add_custom_target(${subdir}_headercheck)
    foreach(header ${headerlist})
      list(FIND dxt_ignore_header "${header}" _index)
      if(${_index} GREATER -1)
        continue()
      endif() # do some name conversion
      set(targname ${header})
      dxt_headercheck_target_name(${targname})
      list(APPEND ${subdir}_dxt_headercheck_targets "${targname}")
      add_dependencies(${subdir}_headercheck ${targname})
    endforeach(header ${headerlist})
  endif(ENABLE_HEADERCHECK)
endmacro(GET_HEADERCHECK_TARGETS)

# cmake-lint: disable=R0915
macro(ADD_SUBDIR_TESTS subdir)
  set(link_xt_libs dunext)
  list(APPEND dxt_test_dirs ${subdir})
  file(GLOB_RECURSE test_sources "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/*.cc")
  foreach(source ${test_sources})
    set(ranks "1")
    if(source MATCHES "mpi")
      list(APPEND ranks ${DUNE_MAX_TEST_CORES})
    endif(source MATCHES "mpi")
    get_filename_component(testbase ${source} NAME_WE)
    string(REPLACE ".cc" ".mini" minifile ${source})
    if(EXISTS ${minifile})
      if(dune-testtools_FOUND)
        dune_add_system_test(
          SOURCE
          ${source}
          ${COMMON_HEADER}
          INIFILE
          ${minifile}
          BASENAME
          test_${testbase}
          CREATED_TARGETS
          targetlist_${testbase}
          ADDED_TESTS
          testlist_${testbase}
          SCRIPT
          dune_xt_execute.py
          ${DEBUG_MACRO_TESTS})
        foreach(target ${targetlist_${testbase}})
          target_link_libraries(${target} ${link_xt_libs} ${COMMON_LIBS} ${GRID_LIBS} gtest_dune_xt)
          list(APPEND ${subdir}_dxt_test_binaries ${target})
          set(dxt_test_names_${target} ${testlist_${testbase}_${target}})
          foreach(test_name ${dxt_test_names_${target}})
            set_tests_properties(${test_name} PROPERTIES LABELS ${subdir})
          endforeach()
        endforeach(target)
      else(dune-testtools_FOUND)
        message("-- missing dune-testtools, disabling test ${source}")
      endif(dune-testtools_FOUND)
    else(EXISTS ${minifile})
      set(target test_${testbase})
      dune_add_test(
        NAME
        ${target}
        SOURCES
        ${source}
        ${COMMON_HEADER}
        LINK_LIBRARIES
        ${link_xt_libs}
        ${COMMON_LIBS}
        ${GRID_LIBS}
        gtest_dune_xt
        COMMAND
        ${CMAKE_BINARY_DIR}/run-in-dune-env
        CMD_ARGS
        ${CMAKE_CURRENT_BINARY_DIR}/${target}
        --gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/${target}.xml
        TIMEOUT
        ${DXT_TEST_TIMEOUT}
        MPI_RANKS
        ${ranks})
      list(APPEND ${subdir}_dxt_test_binaries ${target})
      set(dxt_test_names_${target} ${target})
      set_tests_properties(${target} PROPERTIES LABELS ${subdir})
    endif(EXISTS ${minifile})
  endforeach(source)
  file(GLOB_RECURSE test_templates "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/*.tpl")
  foreach(template ${test_templates})
    set(ranks "1")
    if(template MATCHES "mpi")
      list(APPEND ranks ${DUNE_MAX_TEST_CORES})
    endif(template MATCHES "mpi")
    get_filename_component(testbase ${template} NAME_WE)
    string(REPLACE ".tpl" ".py" config_fn "${template}")
    string(REPLACE ".tpl" ".tpl.cc" out_fn "${template}")
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" out_fn "${out_fn}")
    # get the last completed cache for the codegen execution during configure time
    foreach(mod ${ALL_DEPENDENCIES})
      dune_module_path(MODULE ${mod} RESULT ${mod}_binary_dir BUILD_DIR)
      if(IS_DIRECTORY ${${mod}_binary_dir})
        set(last_dep_bindir ${${mod}_binary_dir})
      endif()
    endforeach(mod DEPENDENCIES)

    dune_execute_process(
      COMMAND
      ${CMAKE_BINARY_DIR}/run-in-dune-env
      dxt_code_generation.py
      "${config_fn}"
      "${template}"
      "${CMAKE_BINARY_DIR}"
      "${out_fn}"
      "${last_dep_bindir}"
      OUTPUT_VARIABLE
      codegen_output)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/codegen.${testbase}.log" ${codegen_output})
    file(GLOB generated_sources "${out_fn}.*")
    if("" STREQUAL "${generated_sources}")
      set(generated_sources ${out_fn})
    endif()
    add_custom_command(
      OUTPUT "${generated_sources}"
      COMMAND ${CMAKE_BINARY_DIR}/run-in-dune-env dxt_code_generation.py "${config_fn}" "${template}"
              "${CMAKE_BINARY_DIR}" "${out_fn}" "${last_dep_bindir}"
      DEPENDS "${config_fn}" "${template}"
      VERBATIM USES_TERMINAL)
    foreach(gen_source ${generated_sources})
      string(REPLACE "${out_fn}." "" postfix "${gen_source}")
      string(REPLACE "${out_fn}" "" postfix "${postfix}")
      string(REPLACE ".cc" "" postfix "${postfix}")
      if(NOT "" STREQUAL "${postfix}")
        set(postfix "__${postfix}")
      endif()
      set(target test_${testbase}${postfix})
      dune_add_test(
        NAME
        ${target}
        SOURCES
        ${gen_source}
        ${COMMON_HEADER}
        LINK_LIBRARIES
        ${link_xt_libs}
        ${COMMON_LIBS}
        ${GRID_LIBS}
        gtest_dune_xt
        COMMAND
        ${CMAKE_BINARY_DIR}/run-in-dune-env
        CMD_ARGS
        ${CMAKE_CURRENT_BINARY_DIR}/${target}
        --gtest_output=xml:${CMAKE_CURRENT_BINARY_DIR}/${target}.xml
        TIMEOUT
        ${DXT_TEST_TIMEOUT}
        MPI_RANKS
        ${ranks})
      list(APPEND ${subdir}_dxt_test_binaries ${target})
      set(dxt_test_names_${target} ${target})
      set_tests_properties(${target} PROPERTIES LABELS ${subdir})
    endforeach()
  endforeach(template ${test_templates})
  add_custom_target(${subdir}_test_templates SOURCES ${test_templates})

  # this excludes meta-ini variation test cases because there binary name != test name
  foreach(test ${${subdir}_xt_test_binaries})
    if(TARGET test)
      set_tests_properties(${test} PROPERTIES TIMEOUT ${DXT_TEST_TIMEOUT})
      set_tests_properties(${test} PROPERTIES LABELS ${subdir})
    endif(TARGET test)
  endforeach()

  add_custom_target(${subdir}_test_binaries DEPENDS ${${subdir}_dxt_test_binaries})
  if("${subdir}" STREQUAL "functions")
    # There are a lot of tests in the function subdir and these take a long time in CI. We thus create two additional
    # targets here, each containing half of the tests.
    list(LENGTH functions_dxt_test_binaries num_tests)
    math(EXPR half_num_tests "${num_tests}/2")
    list(SUBLIST functions_dxt_test_binaries 0 ${half_num_tests} functions1_dxt_test_binaries)
    # If length is -1, all list elements starting from <begin> will be returned
    list(SUBLIST functions_dxt_test_binaries ${half_num_tests} -1 functions2_dxt_test_binaries)
    add_custom_target(functions1_test_binaries DEPENDS ${functions1_dxt_test_binaries})
    add_custom_target(functions2_test_binaries DEPENDS ${functions2_dxt_test_binaries})
    foreach(label functions1 functions2)
      foreach(test ${${label}_dxt_test_binaries})
        set_tests_properties(${test} PROPERTIES LABELS ${label})
      endforeach()
    endforeach()
  endif()
  add_custom_target(
    ${subdir}_check
    COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${DXT_TEST_TIMEOUT} -j ${DXT_TEST_PROCS}
    DEPENDS ${subdir}_test_binaries
    USES_TERMINAL)
  add_custom_target(
    ${subdir}_recheck
    COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${DXT_TEST_TIMEOUT} --rerun-failed -j ${DXT_TEST_PROCS}
    DEPENDS ${subdir}_test_binaries
    USES_TERMINAL)
  foreach(target ${${subdir}_dxt_test_binaries})
    set(all_sorted_testnames "${all_sorted_testnames}/${dxt_test_names_${target}}")
  endforeach()
  set(${subdir}_dxt_headercheck_targets "")
  get_headercheck_targets(${subdir})
endmacro(ADD_SUBDIR_TESTS)

macro(FINALIZE_TEST_SETUP)
  set(combine_targets test_templates test_binaries check recheck)
  foreach(target ${combine_targets})
    add_custom_target(${target})
    foreach(subdir ${dxt_test_dirs})
      add_dependencies(${target} ${subdir}_${target})
    endforeach()
  endforeach()

  foreach(subdir ${dxt_test_dirs})
    set(dxt_test_binaries "${dxt_test_binaries} ${${subdir}_dxt_test_binaries}")
  endforeach()
  # set(${subdir}_dxt_headercheck_targets "")

  if(ALBERTA_FOUND)
    foreach(test ${dxt_test_binaries})
      if(${test} MATCHES alberta_1d)
        add_dune_alberta_flags(GRIDDIM 1 ${test})
      elseif(${test} MATCHES alberta_2d)
        add_dune_alberta_flags(GRIDDIM 2 ${test})
      elseif(${test} MATCHES alberta_3d)
        add_dune_alberta_flags(GRIDDIM 3 ${test})
      endif()
    endforeach()

    foreach(test ${dxt_test_binaries})
      if(${test} MATCHES 2d_simplex_alberta)
        add_dune_alberta_flags(GRIDDIM 2 ${test})
      elseif(${test} MATCHES 3d_simplex_alberta)
        add_dune_alberta_flags(GRIDDIM 3 ${test})
      endif()
    endforeach()
  endif()
endmacro()

macro(DXT_EXCLUDE_FROM_HEADERCHECK)
  exclude_from_headercheck(${ARGV0}) # make this robust to argument being passed with or without ""
  string(REGEX REPLACE "[\ \n]+([^\ ])" ";\\1" list ${ARGV0})
  set(list "${list};${ARGV}")
  foreach(item ${list})
    set(item ${CMAKE_CURRENT_SOURCE_DIR}/${item})
    list(APPEND dxt_ignore_header ${item})
  endforeach()
endmacro(DXT_EXCLUDE_FROM_HEADERCHECK)

macro(DXT_ADD_PYTHON_TESTS)
  add_custom_target(
    xt_test_python
    "${CMAKE_BINARY_DIR}/run-in-dune-env" "py.test" "${CMAKE_BINARY_DIR}/python" "--cov" "${CMAKE_CURRENT_SOURCE_DIR}/"
    "--junitxml=pytest_results.xml"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    DEPENDS bindings
    VERBATIM USES_TERMINAL)
  if(NOT TARGET test_python)
    add_custom_target(test_python)
  endif(TARGET test_python)
  add_dependencies(test_python xt_test_python)
endmacro(DXT_ADD_PYTHON_TESTS)
