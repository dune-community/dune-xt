// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   René Fritze     (2015 - 2016, 2018 - 2019)

#ifndef DUNE_XT_COMMON_LOCALIZATION_STUDY_HH
#define DUNE_XT_COMMON_LOCALIZATION_STUDY_HH

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <dune/common/dynvector.hh>

#include <dune/xt/common/exceptions.hh>

namespace Dune {
namespace XT {
namespace Common {

class LocalizationStudy
{
public:
  LocalizationStudy(const std::vector<std::string> only_these_indicators = {});

  virtual ~LocalizationStudy();

  virtual std::string identifier() const = 0;

  virtual DynamicVector<double> compute_reference_indicators() const = 0;

  virtual std::vector<std::string> provided_indicators() const = 0;

  virtual DynamicVector<double> compute_indicators(const std::string type) const = 0;

  std::vector<std::string> used_indicators() const;

  /*std::map< std::string, std::vector< double > >*/ void run(std::ostream& out = std::cout);

private:
  const std::vector<std::string> only_these_indicators_;
}; // class LocalizationStudy

} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_LOCALIZATION_STUDY_HH
