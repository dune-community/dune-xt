# This script assumes that clang-tidy is installed and in the PATH.
# In addition, clang-tidy needs a compilation database, so you have to run
# cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .
# in your build directory and then copy or symlink the resulting
# compile_commands.json to the source directory.

# It is possible to run clang-tidy on the whole directory at once, using parallelism and all.
# However, even with the run-clang-tidy.py script and -j1, the fixes are applied several
# times to each header (once for each compilation unit in which the header is included),
# resulting in messed up files. We thus simply run clang-tidy sequentially for each
# compilation unit. This takes a long time, but since we only have to do this every
# now and then and can do other stuff while waiting for the script, this should be fine.

# Headerchecks
BUILD_DIR=../build/clang-debug/dune-xt
CLANG_TIDY_BINARY=clang-tidy
for file in `find ${BUILD_DIR}/headercheck -name *.hh.cc`; do
    ${CLANG_TIDY_BINARY} --config-file=.clang-tidy -fix ${file}
done
# Regular C++ files
for file in `find dune/xt/ -name *.cc`; do
    ${CLANG_TIDY_BINARY} --config-file=.clang-tidy -fix ${file}
done
# Python bindings
for file in `find python/ -name *.cc`; do
    ${CLANG_TIDY_BINARY} --config-file=.clang-tidy -fix ${file}
done
