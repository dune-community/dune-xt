// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_LA_SOLVER_DUNEDYNAMIC_HH
#define DUNE_STUFF_LA_SOLVER_DUNEDYNAMIC_HH

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configtree.hh>

#include <dune/stuff/la/container/dunedynamic.hh>

#include "../solver.hh"

namespace Dune {
namespace Stuff {
namespace LA {


template <class S>
class Solver<DuneDynamicMatrix<S>> : protected SolverUtils
{
public:
  typedef DuneDynamicMatrix<S> MatrixType;

  Solver(const MatrixType& matrix)
    : matrix_(matrix)
  {
  }

  static std::vector<std::string> options()
  {
    return {"superlu"};
  }

  static Common::ConfigTree options(const std::string& type)
  {
    SolverUtils::check_given(type, options());
    return Common::ConfigTree("type", type);
  } // ... options(...)

  size_t apply(const DuneDynamicVector<S>& rhs, DuneDynamicVector<S>& solution) const
  {
    return apply(rhs, solution, options()[0]);
  }

  size_t apply(const DuneDynamicVector<S>& rhs, DuneDynamicVector<S>& solution, const std::string& type) const
  {
    this->check_given(type, options());
    // solve
    try {
      matrix_.backend().solve(solution.backend(), rhs.backend());
    } catch (FMatrixError&) {
      return 1;
    }
// check
#ifndef NDEBUG
    if (!rhs.almost_equal(matrix_.mv(solution)))
      return 4;
#endif // NDEBUG
    return 0;
  } // ... apply(...)

  size_t apply(const DuneDynamicVector<S>& rhs, DuneDynamicVector<S>& solution, const Common::ConfigTree& opts) const
  {
    if (!opts.has_key("type"))
      DUNE_THROW_COLORFULLY(Exception::configuration_error,
                            "Given options (see below) need to have at least the key 'type' set!\n\n" << opts);
    const auto type = opts.get<std::string>("type");
    return apply(rhs, solution, type);
  } // ... apply(...)

private:
  const MatrixType& matrix_;
}; // class Solver< DuneDynamicMatrix< ... > >


} // namespace LA
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_LA_SOLVER_DUNEDYNAMIC_HH
