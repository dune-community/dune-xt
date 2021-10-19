macro(dune_xt_module_version_from_git)

    execute_process(COMMAND
        ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/cmake/modules/versioneer.py
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

    message(STATUS "GIT_CWD ${PROJECT_SOURCE_DIR} -- ${GIT_DESCRIBE_VERSION}")

    set(DUNE_MOD_VERSION ${GIT_DESCRIBE_VERSION})
    # Reset variables from dune-common/cmake/modules/DuneMacros.cmake:dune_module_information
    extract_major_minor_version("${DUNE_MOD_VERSION}" DUNE_VERSION)
    set(${DUNE_MOD_NAME_UPPERCASE}_VERSION          "${DUNE_MOD_VERSION}")
    set(${DUNE_MOD_NAME_UPPERCASE}_VERSION_MAJOR    "${DUNE_VERSION_MAJOR}")
    set(${DUNE_MOD_NAME_UPPERCASE}_VERSION_MINOR    "${DUNE_VERSION_MINOR}")
    set(${DUNE_MOD_NAME_UPPERCASE}_VERSION_REVISION "${DUNE_VERSION_REVISION}")

    # configure CPack
    set(CPACK_PACKAGE_NAME "${DUNE_MOD_NAME}")
    set(CPACK_PACKAGE_VERSION "${DUNE_VERSION_MAJOR}.${DUNE_VERSION_MINOR}.${DUNE_VERSION_REVISION}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}" "\\\\.svn" "\\\\.git" ".*/*\\\\.gitignore")

endmacro(dune_xt_module_version_from_git)
