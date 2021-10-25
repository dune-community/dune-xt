macro(dune_xt_module_version_from_git TARGET_MODULE)

    if (dune-xt_MODULE_PATH)
      set(VERSIONEER_DIR ${dune-xt_MODULE_PATH})
    else()
      set(VERSIONEER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
    endif()
    execute_process(COMMAND
        ${PYTHON_EXECUTABLE}
            ${VERSIONEER_DIR}/versioneer.py
            ${PROJECT_SOURCE_DIR}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
                    ERROR_VARIABLE GIT_DESCRIBE_ERROR
                    RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
    if(GIT_DESCRIBE_ERROR_CODE)
      message(FATAL_ERROR "Extracting version information failed: ${GIT_DESCRIBE_ERROR}")
    endif()

    set(${TARGET_MODULE}_VERSION ${GIT_DESCRIBE_VERSION})
    # Reset variables from dune-common/cmake/modules/DuneMacros.cmake:dune_module_information
    extract_major_minor_version("${GIT_DESCRIBE_VERSION}" DUNE_VERSION)
    set(${TARGET_MODULE}_VERSION_MAJOR    "${DUNE_VERSION_MAJOR}")
    set(${TARGET_MODULE}_VERSION_MINOR    "${DUNE_VERSION_MINOR}")
    set(${TARGET_MODULE}_VERSION_REVISION "${DUNE_VERSION_REVISION}")

    set(CPACK_PACKAGE_NAME "${DUNE_MOD_NAME}")
    set(CPACK_PACKAGE_VERSION "${DUNE_VERSION_MAJOR}.${DUNE_VERSION_MINOR}.${DUNE_VERSION_REVISION}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "\\\\.svn" "\\\\.git" ".*/*\\\\.gitignore")

    # reset variables from dune-common/cmake/modules/DuneMacros.cmake:dune_project
    set(ProjectVersion         "${GIT_DESCRIBE_VERSION}")
    set(ProjectVersionString   "${DUNE_VERSION_MAJOR}.${DUNE_VERSION_MINOR}.${DUNE_VERSION_REVISION}")
    set(ProjectVersionMajor    "${DUNE_VERSION_MAJOR}")
    set(ProjectVersionMinor    "${DUNE_VERSION_MINOR}")
    set(ProjectVersionRevision "${DUNE_VERSION_REVISION}")

    # need to re-run template insertion
    configure_file(
      ${CONFIG_VERSION_FILE}
      ${PROJECT_BINARY_DIR}/${ProjectName}-config-version.cmake @ONLY)

endmacro(dune_xt_module_version_from_git)
