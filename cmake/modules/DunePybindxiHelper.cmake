function(dune_pybindxi_add_helper_lib target_name)
  add_library(${target_name} EXCLUDE_FROM_ALL ${ARGN})
  target_include_directories(${target_name}
    PRIVATE ${PYBIND11_INCLUDE_DIR}
    PRIVATE ${PYTHON_INCLUDE_DIRS})
endfunction()
