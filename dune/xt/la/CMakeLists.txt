set(lib_dune_xt_la_sources
    pattern.cc)

dune_library_add_sources(dunextla SOURCES ${lib_dune_la_common_sources})

add_subdirectory(test EXCLUDE_FROM_ALL)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
