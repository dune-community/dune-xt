// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2012 - 2015)
//   Tobias Leibner  (2014)

#include <dune/xt/test/main.hxx>

#include <fstream>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/grid/output/pgf.hh>

using namespace Dune::XT::Grid;

/** output files are compiled in test-compile-pgfoutput.sh target
 **/
TEST(PgfOutput, Sgrid)
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
