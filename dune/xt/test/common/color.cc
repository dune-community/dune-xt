// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014 - 2017)
//   René Fritze     (2012 - 2016, 2018)
//   Stefan Girke    (2012)
//   Tobias Leibner  (2014, 2016, 2018)

#include <dune/xt/common/test/main.hxx>
#include <dune/xt/common/color.hh>

using namespace Dune::XT::Common;

GTEST_TEST(Color, All)
{
  std::cout << "Testing some color codes for this terminal." << std::endl;
  std::cout << "WARNING: This test will succeed although some of the color codes" << std::endl;
  std::cout << "are not supported by your terminal!" << std::endl;
  std::cout << "(Sometimes 'blink', 'reverse' or 'italic' are not supported.)" << std::endl;
  std::cout << StreamModifiers::underline << "a simple 'hello world':" << StreamModifiers::endunderline << " "
            << "hello world" << std::endl;
  std::cout << StreamModifiers::underline << "a colored 'hello world':" << StreamModifiers::endunderline << " "
            << highlight_string("hello world", 1) << std::endl;
  std::cout << StreamModifiers::underline << "a colored 'hello world':" << StreamModifiers::endunderline << " "
            << Colors::bgreen << "hello world" << StreamModifiers::normal << std::endl;
  std::cout << StreamModifiers::underline << "a blinking 'hello world':" << StreamModifiers::endunderline << " "
            << StreamModifiers::blink << "hello world" << StreamModifiers::endblink << std::endl;
  std::cout << StreamModifiers::underline << "an italic 'hello world':" << StreamModifiers::endunderline << " "
            << StreamModifiers::italic << "hello world" << StreamModifiers::enditalic << std::endl;
  std::cout << StreamModifiers::underline << "an underlined 'hello world':" << StreamModifiers::endunderline << " "
            << StreamModifiers::underline << "hello world" << StreamModifiers::endunderline << std::endl;
  std::cout << StreamModifiers::underline << "a reverse 'hello world':" << StreamModifiers::endunderline << " "
            << StreamModifiers::reverse << "hello world" << StreamModifiers::endreverse << std::endl;
  std::cout << StreamModifiers::underline
            << "a 'hello world' with highlighted substrings ('o'):" << StreamModifiers::endunderline << " "
            << highlight_search_string("hello world", "o", 3) << std::endl;
  std::cout << StreamModifiers::underline << "a highlighted 'hello world'-template:" << StreamModifiers::endunderline
            << " " << highlight_template("Hello< World, Hello< World, Hello< World< Hello, World > > > >") << std::endl;
  std::cout << StreamModifiers::underline
            << "a highlighted 'hello world'-template only showing two levels:" << StreamModifiers::endunderline << " "
            << highlight_template("Hello< World, Hello< World, Hello< World< Hello, World > > > >", 2) << std::endl;
  std::cout << StreamModifiers::underline
            << "colored 'hello world' for all available colors( 0 - 255):" << StreamModifiers::endunderline << " "
            << std::endl;
  for (size_t i = 0; i < 256; ++i)
    std::cout << highlight_string("hello world - ", i);
  std::cout << std::endl;
} // main
