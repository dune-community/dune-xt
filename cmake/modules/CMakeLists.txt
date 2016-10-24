set(modules
    DunePybindxiMacros.cmake
    DunePybindxiUtils.cmake
    FindEigen3.cmake
    FindPythonLibsNew.cmake)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
