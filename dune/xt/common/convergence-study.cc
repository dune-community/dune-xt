// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2019)
//   René Fritze     (2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2019 - 2020)

#include "config.h"

#include <algorithm>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/color.hh>

#include "convergence-study.hh"

namespace Dune::XT::Common {


double ConvergenceStudy::expected_rate(const std::string& /*type*/, const std::string& /*id*/) const
{
  return 1;
}

std::vector<std::string> ConvergenceStudy::filter(const std::vector<std::string>& vec,
                                                  const std::vector<std::string>& only_these) const
{
  if (only_these.empty())
    return vec;
  std::vector<std::string> ret;
  for (const auto& element : vec) {
    if (std::find(only_these.begin(), only_these.end(), element) != only_these.end())
      ret.emplace_back(element);
  }
  return ret;

} // ... filter(...)

std::vector<std::pair<std::string, std::string>>
ConvergenceStudy::filter(const std::vector<std::pair<std::string, std::string>>& vec,
                         const std::vector<std::string>& only_these) const
{
  if (only_these.empty())
    return vec;
  std::vector<std::pair<std::string, std::string>> ret;
  for (const auto& element : vec) {
    const auto& key = element.first;
    const auto& value = element.second;
    if (std::find(only_these.begin(), only_these.end(), key) != only_these.end())
      ret.emplace_back(std::make_pair(key, value));
  }
  return ret;

} // ... filter(...)

std::string ConvergenceStudy::lfill(const std::string& id, const size_t len) const
{
  if (id.size() == len)
    return id;
  if (id.size() > len)
    return id.substr(0, len - 1) + ".";
  else /*if (id.size() < len_)*/
    return std::string(len - id.size(), ' ') + id;
}

std::string ConvergenceStudy::cfill(const std::string& id, const size_t len) const
{
  if (id.size() > len)
    return id.substr(0, len - 1) + ".";
  const size_t rpadd = std::ceil((len - id.size()) / 2);
  return std::string(len - rpadd - id.size(), ' ') + id + std::string(rpadd, ' ');
}

double ConvergenceStudy::extract(const std::map<std::string, std::map<std::string, double>>& level_data,
                                 const std::string& type,
                                 const std::string& id) const
{
  try {
    const auto& type_data = level_data.at(type);
    try {
      return type_data.at(id);
    } catch (const std::out_of_range&) {
      DUNE_THROW(Exceptions::requirements_not_met, "data missing for\n   type = " << type << "\n   id = " << id);
    }
  } catch (const std::out_of_range&) {
    DUNE_THROW(Exceptions::requirements_not_met, "data missing for\n   type = " << type);
  }
} // ... extract(...)

double ConvergenceStudy::extract(const std::map<size_t, std::map<std::string, std::map<std::string, double>>>& data,
                                 const size_t level,
                                 const std::string& type,
                                 const std::string& id) const
{
  try {
    return extract(data.at(level), type, id);
  } catch (const std::out_of_range&) {
    DUNE_THROW(Exceptions::requirements_not_met, "data missing for level " << level << "!");
  }
} // ... extract(...)

void ConvergenceStudy::print_eoc(std::ostream& out,
                                 const size_t len,
                                 const std::map<size_t, std::map<std::string, std::map<std::string, double>>>& data,
                                 const size_t level,
                                 const std::string& type,
                                 const std::string& id,
                                 const std::string& target_id) const
{
  auto& self = *this;
  const double quantity_old = extract(data, level - 1, type, id);
  if (FloatCmp::eq(quantity_old, 0.))
    out << lfill("inf", len);
  else {
    const double quantity_new = extract(data, level, type, id);
    const auto target_old = extract(data, level - 1, "target", target_id);
    const auto target_new = extract(data, level, "target", target_id);
    const double eoc_value = std::log(quantity_new / quantity_old) / std::log(target_new / target_old);
    std::stringstream eoc_str;
    if (eoc_value < -999)
      eoc_str << "-inf";
    else if (eoc_value > 9999)
      eoc_str << "inf";
    else
      eoc_str << std::setprecision(static_cast<int>(len - /*dot*/ 1 - /*sign*/ (eoc_value > 0 ? 0 : 1)
                                                    - /*prefix*/ std::ceil(std::abs(std::log10(std::abs(eoc_value))))))
              << std::fixed << eoc_value;
    // color string
    if (eoc_value > (0.9 * self.expected_rate(type, id)))
      out << color_string(lfill(eoc_str.str(), len), Colors::green);
    else if (eoc_value > 0.0)
      out << color_string(lfill(eoc_str.str(), len), Colors::brown);
    else
      out << color_string(lfill(eoc_str.str(), len), Colors::red);
  }
} // ... print_eoc(...)

std::map<std::string, std::map<std::string, std::map<size_t, double>>>
ConvergenceStudy::run(const std::vector<std::string>& only_these, std::ostream& out)
{
  auto& self = *this;
  // check what we want to compute and print
  DUNE_THROW_IF(self.targets().empty(),
                Exceptions::you_are_using_this_wrong,
                "there has to be at least a single target!\n   self.targets() = " << self.targets());
  const std::vector<std::string> actual_targets =
      filter(self.targets(), only_these).empty() ? self.targets() : filter(self.targets(), only_these);
  const auto actual_norms = filter(self.norms(), only_these);
  const auto actual_estimates = filter(self.estimates(), only_these);
  const auto actual_quantities = filter(self.quantities(), only_these);
  DUNE_THROW_IF(actual_norms.size() /*+ actual_estimates.size()*/ + actual_quantities.size() == 0,
                Exceptions::you_are_using_this_wrong,
                "there has to be at least a single norm estimate, or quantity!\n   self.norms() = "
                    << self.norms()
                    //                    << "\n   self.estimates() = "
                    //                    << self.estimates()
                    << "\n   self.quantities() = " << self.quantities() << "\n   only_these = " << only_these);
  const size_t column_width = 8;
  const size_t eoc_column_width = (actual_targets.size() > 1) ? 8 : 4;
  // build header
  const auto disc_info_title = self.discretization_info_title();
  // - discretization
  std::string h1 = " ";
  std::string d1(disc_info_title.size() + 2, '-');
  std::string h2 = " " + disc_info_title + " ";
  std::string delim = d1;
  for (const auto& id : actual_targets) {
    h2 += "| " + lfill(id, column_width) + " ";
    d1 += std::string(column_width + 3, '-');
    delim += "+" + std::string(column_width + 2, '-');
  }
  h1 += cfill("discretization", disc_info_title.size() + actual_targets.size() * (column_width + 3)) + " ";
  // - norms
  for (const auto& id : actual_norms) {
    h1 += "| " + cfill(id, column_width + actual_targets.size() * (eoc_column_width + 3)) + " ";
    d1 += "+" + std::string(column_width + 2 + actual_targets.size() * (eoc_column_width + 3), '-');
    h2 += "| " + lfill("value", column_width) + " ";
    delim += "+" + std::string(column_width + 2, '-');
    for (const auto& target_id : actual_targets) {
      h2 += "| " + lfill(actual_targets.size() == 1 ? "EOC" : "EOC(" + target_id + ")", eoc_column_width) + " ";
      delim += "+" + std::string(eoc_column_width + 2, '-');
    }
  }
// - estimates
#if 0
    for id, _ in actual_estimates:
        if print_full_estimate:
            h1 += "| " + cfill(id, column_width + (len(actual_targets) + 1)*(eoc_column_width + 3)) + " "
        else:
            h1 += "| " + cfill(id, (len(actual_targets) + 1)*(eoc_column_width + 3) - 3) + " "
        if print_full_estimate:
            d1 += "+" + "-"*(column_width + 2)
            h2 += "| " + lfill("estimate", column_width) + " "
            delim += "+" + "-"*(column_width + 2)
        d1 += "+" + "-"*(eoc_column_width + 2)
        h2 += "| " + lfill("eff.", eoc_column_width) + " "
        delim += "+" + "-"*(eoc_column_width + 2)
        for target_id in actual_targets:
            d1 += "+" + "-"*(eoc_column_width + 2)
            h2 += "| " + lfill("EOC" if len(actual_targets) == 1 else "EOC ({})".format(target_id), eoc_column_width) + " "
            delim += "+" + "-"*(eoc_column_width + 2)
#endif // 0
  // - quantities
  for (const auto& id : actual_quantities) {
    std::string first_row;
    std::string second_row;
    std::string third_row;
    auto words = tokenize(id, " ");
    if (id.size() <= column_width) {
      first_row = std::string(column_width, ' ');
      second_row = lfill(id, column_width);
      third_row = std::string(column_width, ' ');
    } else {
      // lets see if we can fit all words in these rows somehow nicely
      int last_row = 0;
      for (auto word_it = words.begin(); word_it != words.end(); ++word_it) {
        const auto& word = *word_it;
        if (word.size() > column_width)
          break;
        if (last_row <= 1 && word.size() <= (column_width - first_row.size() - (first_row.empty() ? 0 : 1))) {
          first_row += (first_row.empty() ? "" : " ") + word;
          last_row = 1;
          words.erase(word_it);
        } else if (last_row <= 2 && word.size() <= (column_width - second_row.size() - (second_row.empty() ? 0 : 1))) {
          second_row += (second_row.empty() ? "" : " ") + word;
          last_row = 2;
          words.erase(word_it);
        } else if (last_row <= 3 && word.size() <= (column_width - third_row.size() - (third_row.empty() ? 0 : 1))) {
          third_row += (third_row.empty() ? "" : " ") + word;
          last_row = 3;
          words.erase(word_it);
        }
      }
      if (words.empty()) {
        // this worked, align right
        first_row = lfill(first_row, column_width);
        second_row = lfill(second_row, column_width);
        third_row = lfill(third_row, column_width);
      } else {
        // the above did not work/cover all words, now brute force
        if (id.size() <= column_width) {
          first_row = lfill(id, column_width);
          second_row = std::string(column_width, ' ');
          third_row = std::string(column_width, ' ');
        } else if (id.size() <= 2 * column_width) {
          first_row = id.substr(0, column_width);
          second_row = lfill(id.substr(column_width, column_width), column_width);
          third_row = std::string(column_width, ' ');
        } else {
          first_row = id.substr(0, column_width);
          second_row = id.substr(column_width, column_width);
          third_row = lfill(id.substr(2 * column_width), column_width);
        }
      }
    }
    h1 += "| " + first_row + " ";
    d1 += "+ " + second_row + " ";
    h2 += "| " + third_row + " ";
    delim += "+" + std::string(column_width + 2, '-');
  }
  // print header
  std::string tmp = delim;
  std::replace(tmp.begin(), tmp.end(), '-', '=');
  out << std::string(h1.size(), '=') << "\n" << h1 << "\n" << d1 << "\n" << h2 << "\n" << tmp << std::endl;
  std::map<size_t, std::map<std::string, std::map<std::string, double>>> data;
  // run actual study
  for (size_t level = 0; level <= self.num_refinements(); ++level) {
    // compute some discretization statistics
    const auto disc_info = self.discretization_info(level);
    // and print them
    out << " " << cfill(disc_info, disc_info_title.size()) << " " << std::flush;
    // do the actual computation
    data[level] = compute(level, actual_norms, actual_estimates, actual_quantities);
    // and print the results
    // - targets
    for (const auto& id : actual_targets) {
      std::stringstream ss;
      ss << std::setprecision(2) << std::scientific << extract(data, level, "target", id);
      out << "| " << lfill(ss.str(), column_width) << " " << std::flush;
    }
    // - norms
    for (const auto& id : actual_norms) {
      std::stringstream ss;
      ss << std::setprecision(2) << std::scientific << extract(data, level, "norm", id);
      out << "| " << lfill(ss.str(), column_width) << " " << std::flush;
      for (const auto& target_id : actual_targets) {
        if (level == 0)
          out << "| " << lfill("----", eoc_column_width) << " " << std::flush;
        else {
          out << "| ";
          print_eoc(out, eoc_column_width, data, level, "norm", id, target_id);
          out << " " << std::flush;
        }
      }
    }
// - estimates
#if 0
     if not "estimate" in self.data[level]:
         self.data[level]["estimate"] = {}
     for estimate_id, norm_id in actual_estimates:
         self.data[level]["estimate"][estimate_id] = self.compute_estimate(level, estimate_id)
         if not norm_id in self.data[level]["norm"]:
             self.data[level]["norm"][norm_id] = self.compute_norm(level, norm_id)
         if print_full_estimate:
             print("| " + lfill("{:.2e}".format(self.data[level]["estimate"][estimate_id]), column_width) + " ", end="")
             sys.stdout.flush()
         print("| " +
                 lfill("{:.2f}".format(self.data[level]["norm"][norm_id]/self.data[level]["estimate"][estimate_id]),
                     eoc_column_width) + " ", end="")
         sys.stdout.flush()
         for target_id in actual_targets:
             if level == 0:
                 print("| " + lfill("----", eoc_column_width) + " ", end="")
                 sys.stdout.flush()
             else:
                 print("| " + print_eoc(self.data[level - 1]["estimate"][estimate_id],
                                          self.data[level]["estimate"][estimate_id],
                                          level, target_id) + " ", end="")
                 sys.stdout.flush()
#endif // 0
    // - quantities
    for (const auto& id : actual_quantities) {
      std::stringstream ss;
      ss << std::setprecision(2) << std::scientific << extract(data, level, "quantity", id);
      out << "| " << lfill(ss.str(), column_width) << " " << std::flush;
    }
    // end of line
    out << "\n";
    if (level < self.num_refinements())
      out << delim;
    out << std::endl;
  }
  // convert data
  std::map<std::string, std::map<std::string, std::map<size_t, double>>> return_data;
  for (const auto& level : data)
    for (const auto& category : level.second)
      for (const auto& type : category.second)
        return_data[category.first][type.first][level.first] = type.second;
  return return_data;
} // ... run(...)


} // namespace Dune::XT::Common
