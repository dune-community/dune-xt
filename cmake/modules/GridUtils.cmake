# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013 - 2017)
#   René Fritze     (2014 - 2015, 2018 - 2019)
#   Tobias Leibner  (2020)
# ~~~

set(GRIDTYPE "YASPGRID" CACHE STRING "GRIDTYPE")

set(GRIDDIM "2" CACHE STRING "GRIDDIM")

add_definitions("-DGRIDDIM=${GRIDDIM}")
add_definitions("-DGRIDTYPE=${GRIDTYPE}")
add_definitions("-D${GRIDTYPE}")

if(ALUGRID_FOUND)
  add_definitions("-DENABLE_ALUGRID=1")
endif(ALUGRID_FOUND)
if(ALBERTA_FOUND)
  add_definitions("-DENABLE_ALBERTA=1")
endif(ALBERTA_FOUND)
if(UG_FOUND)
  add_definitions("-DENABLE_UG=1")
endif(UG_FOUND)
