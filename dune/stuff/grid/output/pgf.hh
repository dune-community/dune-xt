// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_STUFF_GRID_OUTPUT_PGF_HH
#define DUNE_STUFF_GRID_OUTPUT_PGF_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/stuff/grid/walk.hh>
#include <dune/stuff/grid/walk_functors.hh>
#include <dune/common/static_assert.hh>

#include <ostream>
#include <cstdio>
#include <limits>
#include <array>

namespace Dune {
namespace Stuff {
namespace Grid {

typedef std::array<std::string, 7> TexColorArrayType;
namespace {
const TexColorArrayType texcolors_ = {{"black", "red", "blue", "green", "yellow", "cyan", "magenta"}};
}


//! A means to ensure we have exactly 2 coords for tikz to draw, even in 1D
struct PgfCoordWrapper : FieldVector<double, 2>
{
  template <int wdim>
  PgfCoordWrapper(const FieldVector<double, wdim>& vector)
    : FieldVector<double, 2>(0.0)
  {
    for (size_t i = 0; i < size_t(std::min(wdim, 2)); ++i)
      (*this)[i] = vector[i];
  }
};

/** \brief print a tex(tikz) representation of an entity to stream
 *  \ingroup Latex
 *  \param output the ostream used for output
 *
 *  A functor to be used in a \ref GridWalk
 **/
class PgfEntityFunctor
{
public:
  PgfEntityFunctor(std::ostream& output)
    : file_(output)
  {
  }

  template <class Entity>
  void operator()(const Entity& ent, const int ent_idx)
  {
    const auto& geo = ent.geometry();
    for (std::size_t i = 0; i < geo.corners(); ++i) {
      const PgfCoordWrapper corner(geo.corner(i));
      char buffer[150] = {'\0'};
      const int c =
          std::snprintf(buffer, 150, "\\coordinate(C_%d_%lu) at (%f,%f);\n", ent_idx, i, corner[0], corner[1]);
      assert(c > 0);
      file_ << buffer;
    }
    // \draw (A)--(B)--(C)--cycle;
    file_ << "\\draw ";
    for (std::size_t i = 0; i < geo.corners(); ++i) {
      char buffer[50] = {'\0'};
      std::snprintf(buffer, 50, "(C_%d_%lu)--", ent_idx, i);
      file_ << buffer;
    }
    file_ << "cycle;\n ";
  }

private:
  std::ostream& file_;
};

/** \brief print a tex(tikz) representation of a intersection to stream
 *  \ingroup Latex
 *  \param file the ostream used for output
 *  \param color a latex color definiton used for lines
 *  \param print_entitiy_num whether to draw the Entity's Id number
 *
 *  A functor to be used in a \ref GridWalk
 *  \see Pgf
 **/
class PgfEntityFunctorIntersections
{
public:
  PgfEntityFunctorIntersections(std::ostream& file, const std::string color = "black",
                                const bool printEntityIndex = false)
    : file_(file)
    , color_(color)
    , printEntityIndex_(printEntityIndex)
  {
  }

  template <class EntityType>
  void operator()(const EntityType& ent, int ent_idx)
  {
    maybePrintEntityIndex(ent, ent_idx);
  }

  template <class EntityType, class IntersectionType>
  void operator()(const EntityType& /*ent*/, const IntersectionType& intersection)
  {
    typedef typename IntersectionType::Geometry IntersectionGeometry;
    //    typedef Dune::FieldVector< typename IntersectionGeometry::ctype, IntersectionGeometry::coorddimension>
    //        CoordType;
    PgfCoordWrapper a(intersection.geometry().corner(0));
    PgfCoordWrapper b(intersection.geometry().corner(1));
    char buffer[250] = {'\0'};
    std::snprintf(buffer,
                  250,
                  "\\draw[draw=%s,line width=\\gridlinewidth pt,line cap=round] (%f,%f)--(%f,%f);\n",
                  color_.c_str(),
                  a[0],
                  a[1],
                  b[0],
                  b[1]);
    file_ << buffer;
    file_.flush();
  }

  template <class Entity>
  void maybePrintEntityIndex(const Entity& entity, const int idx)
  {
    if (!printEntityIndex_)
      return;
    //    typedef typename Entity::Geometry
    //        EntityGeometryType;
    PgfCoordWrapper center(entity.geometry().center());
    char buffer[50] = {'\0'};
    std::snprintf(buffer, 50, "\\node[circle] at (%f,%f) {%d};\n", center[0], center[1], idx);
    file_ << buffer;
  }

protected:
  std::ostream& file_;
  const std::string color_;
  const bool printEntityIndex_;
};

/** \brief print a tex(tikz) representation of each entity to strean
 * basically the same as \ref PgfEntityFunctorIntersections, but shifting the drawn line "inwards"
 *  \ingroup Latex
 *  \param file the ostream used for output
 *  \param color a latex color definiton used for lines
 *  \param print_entitiy_num whether to draw the Entity's Id number
 *
 *  A functor to be used in a \ref GridWalk
 **/
class PgfEntityFunctorIntersectionsWithShift : public PgfEntityFunctorIntersections
{
public:
  PgfEntityFunctorIntersectionsWithShift(std::ostream& file, const std::string color = "black", const int level = 0,
                                         bool printEntityIndex = false)
    : PgfEntityFunctorIntersections(file, color, printEntityIndex)
    , level_(level)
  {
  }

  template <class Entity>
  void operator()(const Entity& ent, const int ent_idx)
  {
    maybePrintEntityIndex(ent, ent_idx);
  }

  template <class EntityType, class IntersectionType>
  void operator()(const EntityType& ent, const IntersectionType& intersection)
  {
    const PgfCoordWrapper center(ent.geometry().center());
    const float fac = 0.16 * level_;
    PgfCoordWrapper a(intersection.geometry().corner(0));
    PgfCoordWrapper b(intersection.geometry().corner(1));
    PgfCoordWrapper a_c = center - a;
    PgfCoordWrapper b_c = center - b;
    a_c *= fac;
    b_c *= fac;
    a += a_c;
    b += b_c;
    char buffer[150] = {'\0'};
    std::snprintf(buffer,
                  150,
                  "\\draw[draw=%s,line width=\\gridlinewidth pt,line cap=round] (%f,%f)--(%f,%f);\n",
                  color_.c_str(),
                  a[0],
                  a[1],
                  b[0],
                  b[1]);
    file_ << buffer;
    file_.flush();
  }

private:
  const int level_;
};

/** \brief Provider of three types of latex drawings for \ref Grid
 *  \tparam GridType a \ref Grid implementation
 **/
template <class GridType>
class PgfOutput
{
public:
  PgfOutput(GridType& grid)
    : grid_(grid)
  {
  }

  //! print a tex representation of any given 2D grid leaf level to file
  void leaf(std::ostream& file, const bool includable = true) const
  {
    if (!includable) {
      file << "\\documentclass{article}\n"
              "\\usepackage{tikz}\n"
              "\\usetikzlibrary{calc,intersections, arrows,shapes.misc,shapes.arrows}\n"
              "\\pagestyle{empty}\n\\newcommand{\\gridplotscale}{14.8}\n"
              "\\newcommand{\\gridcoordscale}{0.8}\n"
              "\\newcommand{\\gridlinewidth}{0.8}\n"
              "\\begin{document}\n"
              "\\begin{tikzpicture}[scale=\\gridplotscale]\n";
    } else
      file << "\\begin{tikzpicture}\n";
    GridWalk<typename GridType::LeafGridView> gridWalk(grid_.leafView());
    PgfEntityFunctorIntersections pgf(file);
    gridWalk(pgf, pgf);

    file << "\\end{tikzpicture}\n";
    if (!includable)
      file << "\\end{document}\n";
  }

  /** \brief refine given grid refineLevel times and draw all entities into the same figure
   * \param includable true: only output tikzpicture environment
   *                   false: wrap tikzpicture environment in complete document
   * \param refineLevel passed to \ref Grid::globalRefine
   * \param file the output stream
   **/
  void stacked(std::ostream& file, const int refineLevel, const bool includable = true)
  {
    if (!includable) {
      file << "\\documentclass{article}\n"
              "\\usepackage{tikz}\n"
              "\\usetikzlibrary{calc,intersections, arrows,shapes.misc,shapes.arrows}\n"
              "\\pagestyle{empty}\n\\newcommand{\\gridplotscale}{14.8}\n"
              "\\newcommand{\\gridcoordscale}{0.8}\n"
              "\\newcommand{\\gridlinewidth}{0.8}\n"
              "\\begin{document}\n"
              "\\begin{tikzpicture}[scale=\\gridplotscale]\n";
    } else
      file << "\\begin{tikzpicture}\n";
    grid_.globalRefine(refineLevel);
    for (int i = 0; i < refineLevel; ++i) {
      typedef typename GridType::LevelGridView ViewType;
      const ViewType& view = grid_.levelView(i);
      GridWalk<ViewType> gridWalk(view);
      PgfEntityFunctorIntersectionsWithShift pgf(file, texcolors_[std::min(i, int(texcolors_.size()))], i, true);
      gridWalk(pgf);
      file << "%%%%%%%%%%%%%%%" << view.size(0) << "%%%%%%%%%%%%%%%%\n";
    }

    file << "\\end{tikzpicture}\n";
    if (!includable)
      file << "\\end{document}\n";
  }

  /** \brief refine given grid refineLevel times and create a subfigure for each level
   *  \param includable true: only output figure environment
   *                   false: wrap figure environment in complete document
   *  \param refineLevel passed to \ref Grid::globalRefine
   *  \param file the output stream
   *
   *  The written latex code allows tweaking by changing values of gridplotscale,
   *  gridcoordscale and gridlinewidth
   **/
  void refineseries(std::ostream& file, const int refineLevel, const bool includable = true)
  {
    if (!includable) {
      file << "\\documentclass{article}\n"
              "\\usepackage{tikz}\n\\usepackage{subfig}\n"
              "\\usetikzlibrary{calc,intersections, arrows,shapes.misc,shapes.arrows}\n"
              "\\pagestyle{empty}\n\\newcommand{\\gridplotscale}{0.8}\n"
              "\\newcommand{\\gridcoordscale}{0.8}\n"
              "\\newcommand{\\gridlinewidth}{0.8}\n"
              "\\begin{document}\n";
    }
    file << "\\begin{figure}\n";
    grid_.globalRefine(refineLevel);
    for (int i = 0; i < refineLevel; ++i) {
      typedef typename GridType::LevelGridView ViewType;
      {
        const ViewType& view = grid_.levelView(i);
        char buffer[80] = {'\0'};
        std::snprintf(buffer, 80, "\\subfloat[Level %d]{\n\\begin{tikzpicture}[scale=\\gridplotscale]\n", i);
        file << buffer;
        GridWalk<ViewType> gridWalk(view);
        PgfEntityFunctorIntersections thisLevel(file, "black", true);
        gridWalk(thisLevel, thisLevel);
      }

      GridWalk<typename GridType::LeafGridView> leafWalk(grid_.leafView());
      typedef typename GridType::LeafGridView::Traits::template Codim<0>::Entity EntityType;
      MinMaxCoordinateFunctor<EntityType> minMaxCoord;
      leafWalk(minMaxCoord);

      switch (GridType::dimensionworld) {
        case 1: {
          char buffer[80]     = {'\0'};
          const double offset = 0.2;
          const char* format = "\\node[scale=\\gridcoordscale] at (%f,0) {(%d)};\n";
          std::snprintf(buffer, 80, format, minMaxCoord.minima_[0] - offset, minMaxCoord.minima_[0]);
          file << buffer;
          std::snprintf(buffer, 80, format, minMaxCoord.maxima_[0] - offset, minMaxCoord.maxima_[0]);
          file << buffer;
          break;
        }
        case 2:
        case 3: {
          const double offset = 0.2;
          const char* format  = "\\node[scale=\\gridcoordscale] at (%f,%f) {(%d,%d)};\n";
          char buffer[100] = {'\0'};
          std::snprintf(buffer,
                        100,
                        format,
                        minMaxCoord.minima_[0] - offset,
                        minMaxCoord.minima_[1] - offset,
                        minMaxCoord.minima_[0],
                        minMaxCoord.minima_[1]);
          file << buffer;
          std::snprintf(buffer,
                        100,
                        format,
                        minMaxCoord.maxima_[0] - offset,
                        minMaxCoord.maxima_[1] - offset,
                        minMaxCoord.maxima_[0],
                        minMaxCoord.maxima_[1]);
          file << buffer;
          break;
        }
      }
      if ((i + 1) % 3 == 0) // three subfigs per row
        file << "\\end{tikzpicture}}\\\\\n";
      else
        file << "\\end{tikzpicture}}\n";
    }
    file << "\\end{figure}\n";
    if (!includable)
      file << "\\end{document}\n";
  }

private:
  GridType& grid_;
};

} // namespace Stuff
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_IO_LATEX_PGF_HH

/** Copyright (c) 2012, Rene Milk
   * All rights reserved.
   *
   * Redistribution and use in source and binary forms, with or without
   * modification, are permitted provided that the following conditions are met:
   *
   * 1. Redistributions of source code must retain the above copyright notice, this
   *    list of conditions and the following disclaimer.
   * 2. Redistributions in binary form must reproduce the above copyright notice,
   *    this list of conditions and the following disclaimer in the documentation
   *    and/or other materials provided with the distribution.
   *
   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   * The views and conclusions contained in the software and documentation are those
   * of the authors and should not be interpreted as representing official policies,
   * either expressed or implied, of the FreeBSD Project.
   **/
