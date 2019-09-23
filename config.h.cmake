// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2013 - 2019)
//   Tobias Leibner  (2014, 2016 - 2019)

/* begin dune-xt-common */
// NEVER delete/alter above comment, dune's cmake relies on it

/* Define to the version of dune-xt-common */
#define DUNE_XT_COMMON_VERSION ${DUNE_XT_COMMON_VERSION}

/* Define to the major version of dune-xt-common */
#define DUNE_XT_COMMON_VERSION_MAJOR ${DUNE_XT_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-xt-common */
#define DUNE_XT_COMMON_VERSION_MINOR ${DUNE_XT_COMMON_VERSION_MINOR}

/* Define to the revision of dune-xt-common */
#define DUNE_XT_COMMON_VERSION_REVISION ${DUNE_XT_COMMON_VERSION_REVISION}

#ifndef HAVE_CBLAS
#cmakedefine01 HAVE_CBLAS
#endif

#ifndef HAVE_EIGEN
#define HAVE_EIGEN ENABLE_EIGEN
#endif

#ifndef HAVE_CLP
#cmakedefine01 HAVE_CLP
#endif

#ifndef HAVE_QHULL
#cmakedefine01 HAVE_QHULL
#endif

#ifndef HAVE_LPSOLVE
#cmakedefine01 HAVE_LPSOLVE
#endif

#ifndef HAVE_LAPACKE
#cmakedefine01 HAVE_LAPACKE
#endif

#ifndef HAVE_MKL
#cmakedefine01 HAVE_MKL
#endif

#ifndef HAVE_MATEXP
#cmakedefine01 HAVE_MATEXP
#endif

#ifndef HAVE_TBB
#cmakedefine01 HAVE_TBB
#endif

#ifndef DXT_DISABLE_LARGE_TESTS
#define DXT_DISABLE_LARGE_TESTS 0
#endif

#ifndef HAS_WORKING_UNUSED_ATTRIBUTE
#cmakedefine HAS_WORKING_UNUSED_ATTRIBUTE 1
#endif

#ifndef HAVE_LIKWID
#cmakedefine01 HAVE_LIKWID
#endif

#ifndef ENABLE_PERFMON
#cmakedefine01 ENABLE_PERFMON
#endif

#if ENABLE_PERFMON && HAVE_LIKWID
#define LIKWID_PERFMON 1
#endif

#ifndef HAVE_MAP_EMPLACE
#cmakedefine HAVE_MAP_EMPLACE 1
#endif

#ifndef DS_MAX_MIC_THREADS
#define DS_MAX_MIC_THREADS ${DS_MAX_MIC_THREADS}
#endif

#ifndef HAVE_DUNE_FEM_PARAMETER_REPLACE
#define HAVE_DUNE_FEM_PARAMETER_REPLACE 0
#endif

/* needed in dune/xt/common/profiler.hh */
#ifndef DUNE_XT_COMMON_DO_TIMING
#define DUNE_XT_COMMON_DO_TIMING 0
#endif

#cmakedefine01 DUNE_XT_WITH_PYTHON_BINDINGS

/*** Silence implicitly False evaluation of undefined macro warnings ****/
#ifndef HAVE_DUNE_FEM
#define HAVE_DUNE_FEM 0
#endif

#ifndef HAVE_DUNE_GRID_MULTISCALE
#define HAVE_DUNE_GRID_MULTISCALE 0
#endif

#ifndef HAVE_DUNE_GRID_GLUE
#define HAVE_DUNE_GRID_GLUE 0
#endif

#ifndef HAVE_ALBERTA
#define HAVE_ALBERTA 0
#endif

#ifndef ALBERTA_DIM
#define ALBERTA_DIM 2
#endif

#ifndef HAVE_DUNE_SPGRID
#define HAVE_DUNE_SPGRID 0
#endif

#ifndef HAVE_EMPLACE
#define HAVE_EMPLACE 0
#endif

#ifndef HAVE_FASP
#define HAVE_FASP 0
#endif

#ifndef ENABLE_SUPERLU
#define ENABLE_SUPERLU 0
#endif

#ifndef HAVE_SUPERLU
#define HAVE_SUPERLU 0
#endif

#ifndef ENABLE_UMFPACK
#define ENABLE_UMFPACK 0
#endif

#ifndef HAVE_UMFPACK
#define HAVE_UMFPACK 0
#endif

#ifndef HAVE_UG
#define HAVE_UG 0
#endif

#ifndef ENABLE_BOOST
#define ENABLE_BOOST 0
#endif

#ifndef ENABLE_PETSC
#define ENABLE_PETSC 0
#endif

#ifndef ENABLE_MPI
#define ENABLE_MPI 0
#endif

#ifndef HAVE_MPI
#define HAVE_MPI 0
#endif

#ifndef HAVE_CONSTEXPR
#define HAVE_CONSTEXPR 0
#endif

#ifndef HAVE_PETSC
#define HAVE_PETSC 0
#endif

#ifndef DUNE_FEM_COMPATIBILITY
#define DUNE_FEM_COMPATIBILITY 0
#endif

#ifndef HAVE_MAP_EMPLACE
#define HAVE_MAP_EMPLACE 0
#endif

/*** End: Silence implicitly False evaluation of undefined macro warnings ****/

#include <boost/config.hpp>
#if HAVE_TBB && defined(BOOST_CLANG)
  // Hack to fix compilation with clang as tbb does not detect C++11 feature correctly for clang. Recent versions of TBB
  // allow to set the macro TBB_USE_GLIBCXX_VERSION to the proper version of libstdc++ to fix this issue, see
  // https://www.threadingbuildingblocks.org/docs/help/reference/appendices/known_issues/linux_os.html. For older versions
  // we need the hack below.
  #define TBB_USE_GLIBCXX_VERSION 40902
  #include "tbb/tbb_stddef.h"
  #if TBB_INTERFACE_VERSION < 4400
    #include <tbb/tbb_config.h>
    #undef __TBB_CPP11_RVALUE_REF_PRESENT
    #undef __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT
    #undef __TBB_CPP11_DECLTYPE_PRESENT
    #undef __TBB_CPP11_LAMBDAS_PRESENT
    #define __TBB_CPP11_RVALUE_REF_PRESENT 1
    #define __TBB_CPP11_VARIADIC_TEMPLATES_PRESENT 1
    #define __TBB_CPP11_DECLTYPE_PRESENT 1
    #define __TBB_CPP11_LAMBDAS_PRESENT 1
  #endif // TBB_INTERFACE_VERSION < 4400
#endif // HAVE_TBB

// This is an unfortunate hack, see the header for an explanation.
#include <dune/xt/common/fix-ambiguous-std-math-overloads.hh>

/* end dune-xt-common */
// NEVER delete/alter above comment, dune's cmake relies on it
