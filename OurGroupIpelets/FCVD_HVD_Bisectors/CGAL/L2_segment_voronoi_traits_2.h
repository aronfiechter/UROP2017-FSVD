// created by Aron Fiechter on 2017-07-13.
// This file is a copy of L2_voronoi_traits_2, except that it is for segments
// instead of points.
// Actually, it is more similar to Env_sphere_traits.h in the end, because for
// intersections of segments we needed conics.

// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Envelope_3/include/CGAL/Env_plane_traits_3.h $
// $Id: Env_plane_traits_3.h 51989 2009-09-21 10:55:53Z efif $
//
// Author(s)     : Ophir Setter

#ifndef CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
#define CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Env_plane_traits_3_functions.h>

namespace CGAL{

template <class ConicTraits_2>
class L2_segment_voronoi_traits_2 : public ConicTraits_2 {

public:
  typedef ConicTraits_2                                     Traits_2;

  typedef typename Traits_2::Point_2                        Point_2;
  typedef typename Traits_2::Curve_2                        Curve_2;
  typedef typename Traits_2::X_monotone_curve_2             X_monotone_curve_2;
  typedef typename Traits_2::Multiplicity                   Multiplicity;

  typedef typename Traits_2::Rat_kernel                     Rat_kernel;
  typedef typename Traits_2::Alg_kernel                     Alg_kernel;
  typedef typename Traits_2::Nt_traits                      Nt_traits;

  /* For point-segment bisectors */
  typedef typename Traits_2::Curve_2                        Conic_arc_2;

  typedef typename Rat_kernel::FT                           Rational;
  typedef typename Rat_kernel::Point_2                      Rat_point_2;
  typedef typename Rat_kernel::Segment_2                    Rat_segment_2;
  typedef typename Rat_kernel::Line_2                       Rat_line_2;
  // typedef typename Rat_kernel::Circle_2                     Rat_circle_2;
  // typedef typename Rat_kernel::Point_3                      Rat_point_3;
  typedef typename Rat_kernel::Ray_2                        Rat_ray_2;
  typedef typename Rat_kernel::Direction_2                  Rat_direction_2;

  typedef typename Alg_kernel::FT                           Algebraic;
  typedef typename Alg_kernel::Point_2                      Alg_point_2;
  // typedef typename Alg_kernel::Circle_2                     Alg_circle_2;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Rat_segment_2                  Xy_monotone_surface_3;
  typedef Rat_segment_2                  Surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity>       Intersection_curve;

  /* Returns the squared distance between two points in L2 metric. */
  static Rational sqdistance(const Point_2& p1, const Point_2& p2) {
    Rational diffx = p1.x() - p2.x();
    Rational diffy = p1.y() - p2.y();
    Rational sqdist = diffx*diffx + diffy*diffy;
    return sqdist;
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static Rational sqdistance(const Point_2& p, const Rat_segment_2& s){
    /* find projection of p on supporting line of s */
    Rat_line_2 l = s.supporting_line();
    Point_2 proj = l.projection(p);

    /* if the projection is on s, the distance is d(p,proj) */
    if (s.has_on(proj)) {
      return sqdistance(p, proj);
    }
    /* otherwise, the distance is min(d(p,s1),d(p,s2)), where s1 and s2 are the
    endpoints of the segment s */
    else {
      return CGAL::min(sqdistance(p, s.source()), sqdistance(p, s.target()));
    }
  }

public:

  class Make_xy_monotone_3 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool /* is_lower */,
                                OutputIterator o) const {
      /* the surfaces we are considering are distance functions from line
       * segments and are already xy_monotone because there is only one possible
       * distance value for any point on the plane */
      *o++ = s; // just insert the surface in o, return o one past the end
      return o;
    }
  };

  Make_xy_monotone_3 make_xy_monotone_3_object() const {
    return Make_xy_monotone_3();
  }

  class Construct_projected_boundary_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      /* the surfaces we are considering are distance functions from line
       * segments and are infinite, so they have no projected boundary */
      return o; // the iterator remains empty
    }
  };

  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }

  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::EQUAL;
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::EQUAL;
    }

    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      // should happen only if the points are equal.
      return CGAL::EQUAL;
    }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3();
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::EQUAL;
    }

  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::EQUAL;
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }


  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {
      if (s1 == s2) {
        Conic_arc_2   c6 = Conic_arc_2 (
                1, 0, 0, 0, 1, 0,       // The parabola.
                CGAL::CLOCKWISE,
                Point_2 (-1.73, -3),    // Approximation of the source.
                0, 0, 0, 0, 1, 3,       // The line: y = -3.
                Point_2 (1.41, -2),     // Approximation of the target.
                0, 0, 0, 0, 1, 2        // The line: y = -2.
        );
        return o;
      } else {
        *o++ = CGAL::make_object(Intersection_curve(CGAL::bisector(s1, s2), 1));
        return o;
      }
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }

};

} //namespace CGAL

#endif // CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
