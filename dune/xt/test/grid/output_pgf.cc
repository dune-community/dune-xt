// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2012 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2016)

#include <dune/xt/test/main.hxx>

#include <fstream>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/grid/output/pgf.hh>

using namespace Dune::XT::Grid;

/** output files are compiled in test-compile-pgfoutput.sh target
 **/
GTEST_TEST(PgfOutput, YaspGrid)
{
  const size_t dim = 2;
  typedef Dune::YaspGrid<dim> GridType;
  std::array<int, dim> n;
  Dune::FieldVector<double, dim> h;

  for (size_t i = 0; i < dim; ++i) {
    n[i] = 2;
    h[i] = 1.0;
  }
  GridType grid(h, n);
  PgfOutput<GridType> output(grid);
  const int max_refines = 2;
  const bool includable = false;
  std::ofstream fileB("pgfoutput_refineseries.tex");
  output.refineseries(fileB, max_refines, includable);
  std::ofstream fileC("pgfoutput_stacked.tex");
  output.stacked(fileC, max_refines, includable);
  std::ofstream fileA("pgfoutput_leaf.tex");
  output.leaf(fileA, includable);
}
