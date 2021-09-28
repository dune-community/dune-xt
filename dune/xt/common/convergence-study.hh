// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2014, 2016 - 2019)
//   René Fritze     (2014 - 2016, 2018 - 2020)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_COMMON_CONVERGENCE_STUDY_HH
#define DUNE_XT_COMMON_CONVERGENCE_STUDY_HH

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <dune/xt/common/logging.hh>

namespace Dune::XT::Common {


class ConvergenceStudy
{
public:
  virtual ~ConvergenceStudy() = default;

  /**
   * \brief The number of refinements to carry out.
   *
   * discretiztion_info() and compute() will be called for 0 <= refinement_level <= num_refinements.
   */
  virtual size_t num_refinements() const = 0;

  /**
   * \brief The identifiers of the targets with respect to which the EOCs will be computed.
   * \note  This must not be empty!
   * \note  Target ids are shortened to length 4 in EOC column headers.
   *
   * This is usually something like "h" or "dt".
   */
  virtual std::vector<std::string> targets() const = 0;

  /**
   * \brief Error norms of which the EOCs will be computed.
   */
  virtual std::vector<std::string> norms() const
  {
    return {};
  }

  /**
   * \brief A posteriori error estimates (and corresponding norms) of which the EOCs will be computed.
   */
  virtual std::vector<std::pair<std::string, std::string>> estimates() const
  {
    return {};
  }

  /**
   * \brief Quantities to be displayed.
   */
  virtual std::vector<std::string> quantities() const
  {
    return {};
  }

  /**
   * \brief Header to be printed in the column which displays discretization_info().
   */
  virtual std::string discretization_info_title() const = 0;

  /**
   * \brief Statistics about the discretization, in addition to targets().
   * \note  This is the first called on a new level and can be used to set the internal current refinemen_level.
   *
   * This is usually somthing like the number of DoFs or grid elements
   */
  virtual std::string discretization_info(const size_t refinement_level) = 0;

  /**
   * \brief Carries out the actual computations on this refinement_level.
   *
   * This usually solves the problem and computes various statistics, which are expected to be of the form
\code
for (std::string type : {"target", "norm", "estimate"}) {
  for (std::string id : this->type()) { // <- this is pseudocode
    if (id in only_these) {
      data[type][id] = some_quantity;
    }
  }
}
\endcode
   * and to contain data["info"]["time_to_solution"] (in seconds).
   */
  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level,
          const std::vector<std::string>& actual_norms,
          const std::vector<std::pair<std::string, std::string>>& actual_estimates,
          const std::vector<std::string>& actual_quantities) = 0;

  /**
   * \brief The expected rates of convergence, merely to be used to color the EOCs.
   */
  virtual double expected_rate(const std::string& type, const std::string& id) const;

protected:
  // some helpers
  std::vector<std::string> filter(const std::vector<std::string>& vec,
                                  const std::vector<std::string>& only_these) const;

  std::vector<std::pair<std::string, std::string>> filter(const std::vector<std::pair<std::string, std::string>>& vec,
                                                          const std::vector<std::string>& only_these) const;

  std::string lfill(const std::string& id, const size_t len) const;

  std::string cfill(const std::string& id, const size_t len) const;

  /// \todo improve
  double extract(const std::map<std::string, std::map<std::string, double>>& level_data,
                 const std::string& type,
                 const std::string& id) const;

  double extract(const std::map<size_t, std::map<std::string, std::map<std::string, double>>>& data,
                 const size_t level,
                 const std::string& type,
                 const std::string& id) const;

  void print_eoc(std::ostream& out,
                 const size_t len,
                 const std::map<size_t, std::map<std::string, std::map<std::string, double>>>& data,
                 const size_t level,
                 const std::string& type,
                 const std::string& id,
                 const std::string& target_id) const;

public:
  /**
   * \brief Runs the study and displays a table with all targets, norms, estimates and quantities given by the study
   *        (if only_these is empty) or only those which are contained in only_these (else).
   *
   * \return data with data[foo][bar][level] where 0 <= level <= num_refinements(), "foo" is one of
   *         {"target", "norm", "estimate", "quantity"} and bar is one of targets(), norms(), estimates() or
   *         quantities() (intersected with only_these), respectively.
   **/
  std::map<std::string, std::map<std::string, std::map<size_t, double>>>
  run(const std::vector<std::string>& only_these = {}, std::ostream& out = std::cout);
}; // class ConvergenceStudy


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_CONVERGENCE_STUDY_HH
