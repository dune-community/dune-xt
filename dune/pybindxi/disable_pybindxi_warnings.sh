# specifies all pybind11 headers as system headers to disable warnings
# should also work for clang which supports some GCC macros for compatibility
# including as system folder does not work as we would have to include the whole dune-xt folder using -isystem
sed -i '1 i\#pragma GCC system_header' *.h
sed -i '1 i\#pragma GCC system_header' detail/*.h
sed -i '1 i\#pragma GCC system_header' stl/*.h
