// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_LA_SOLVER_ISTL_PRECONDITIONERS_HH
#define DUNE_XT_LA_SOLVER_ISTL_PRECONDITIONERS_HH

#include <type_traits>
#include <cmath>

#include <dune/istl/preconditioners.hh>

namespace Dune::XT::LA {


template <class O>
class IdentityPreconditioner : public Dune::Preconditioner<typename O::domain_type, typename O::range_type>
{
public:
  //! \brief The domain type of the preconditioner.
  using domain_type = typename O::domain_type;
  //! \brief The range type of the preconditioner.
  using range_type = typename O::range_type;
  //! \brief The field type of the preconditioner.
  using field_type = typename range_type::field_type;
  using InverseOperator = O;

  IdentityPreconditioner(const SolverCategory::Category cat)
    : category_(cat)
  {}

  //! Category of the preconditioner (see SolverCategory::Category)
  SolverCategory::Category category() const final
  {
    return category_;
  }

  void pre(domain_type&, range_type&) final {}

  void apply(domain_type& v, const range_type& d) final
  {
    v = d;
  }

  void post(domain_type&) final {}

private:
  SolverCategory::Category category_;
};


} // namespace Dune::XT::LA

#endif // DUNE_XT_LA_SOLVER_ISTL_PRECONDITIONERS_HH
