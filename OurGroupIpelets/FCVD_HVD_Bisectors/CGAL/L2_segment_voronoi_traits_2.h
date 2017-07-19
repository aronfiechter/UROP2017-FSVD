// created by Aron Fiechter on 2017-07-13.
// This file implements a model of the concept EnvelopeTraits_3.
// It is mostly is a copy of L2_voronoi_traits_2, except that it is for segments
// instead of points.

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
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>

#include <CGAL/Parabola_segment_2.h>

namespace CGAL {

template <class Conic_traits_2, class Kernel_>
class L2_segment_voronoi_traits_2 : public Conic_traits_2 {

public:
  typedef Kernel_                                   Kernel;
  typedef Conic_traits_2                            C_traits_2;
  typedef L2_segment_voronoi_traits_2<C_traits_2, Kernel> Self;

  typedef typename C_traits_2::Point_2              Point_2;
  typedef typename C_traits_2::Curve_2              Curve_2;
  typedef typename C_traits_2::X_monotone_curve_2   X_monotone_curve_2;
  typedef typename C_traits_2::Multiplicity         Multiplicity;

  typedef typename C_traits_2::Rat_kernel           Rat_kernel;
  typedef typename C_traits_2::Alg_kernel           Alg_kernel;
  typedef typename C_traits_2::Nt_traits            Nt_traits;

  typedef typename Rat_kernel::FT                   Rational;
  typedef typename Rat_kernel::Point_2              Rat_point_2;
  typedef typename Rat_kernel::Segment_2            Rat_segment_2;
  typedef typename Rat_kernel::Line_2               Rat_line_2;
  typedef typename Rat_kernel::Ray_2                Rat_ray_2;

  typedef typename Alg_kernel::FT                   Algebraic;
  typedef typename Alg_kernel::Point_2              Alg_point_2;
  typedef typename Alg_kernel::Segment_2            Alg_segment_2;
  typedef typename Alg_kernel::Line_2               Alg_line_2;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Alg_segment_2                             Surface_3;
  typedef Surface_3                                 Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity>    Intersection_curve;

  /* Returns the squared distance between two points in L2 metric. */
  static Algebraic sqdistance(const Alg_point_2& p1, const Alg_point_2& p2) {
    Algebraic diffx = p1.x() - p2.x();
    Algebraic diffy = p1.y() - p2.y();
    Algebraic sqdist = CGAL::square(diffx) + CGAL::square(diffy);
    return sqdist;
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static Algebraic sqdistance(const Alg_point_2& p, const Alg_segment_2& s) {
    /* find projection of p on supporting line of s */
    Alg_line_2 l = s.supporting_line();
    Alg_point_2 proj = l.projection(p);

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

  /* Construct a point in the middle of the curve cv. This function is copied
   * from Env_sphere_traits_3.h */
  static Alg_point_2 construct_middle_point(const X_monotone_curve_2& cv) {
    /* get the x-value of the middle point */
    Alg_kernel k;
    Alg_point_2 mid_x = k.construct_midpoint_2_object()(
      cv.source(),
      cv.target()
    );

    /* if cv is vertical, it is just a segment */
    if (cv.is_vertical()) return Alg_point_2(mid_x);
    /* otherwise take the point with the same x coordinate but on cv */
    else return Alg_point_2(cv.point_at_x(mid_x));
  }

  // /* Converts a parabola segment (a parabolic arc) into a series of segments.
  //  * This is needed because we're using a linear Kernel, and because the CGAL
  //  * Ipelet interface does not allow to draw parabolic arcs (yet). */
  // void par_arc_to_segments(const Parabola_segment_2& pb_arc, const std::list<Segment_2>& segments) {
  //   /* generate points on the parabolic arc */
  //   std::vector<Point_2> points;
  //   pb_arc.generate_points(points);
  //
  //   /* create segments from points */
  //   for (auto i = 0; i < points.size() - 1; ++i) {
  //     segments.push_back(Segment_2(points[i], points[i + 1]));
  //   }
  // }



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



  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {
      /* if the two segments are the same, their distance function is the same,
       * so there is no intersection */
      if (s1 == s2) {
        return o;
      }
      /* otherwise, for now just make a random polycurve between the two
       * segments, to test if this works */
      // TODO fake
      else {

        Curve_2 c1(
          1, 0, 0, 0, -1, 0, CGAL::COUNTERCLOCKWISE,
          Point_2(Algebraic(0), Algebraic(0)),
          Point_2(Algebraic(3), Algebraic(9))
        );
        X_monotone_curve_2 mc1(c1);
        *o++ = CGAL::make_object(
          Intersection_curve(mc1, 0)
        );
        return o;
      }

    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }



  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      Alg_point_2 ap(p);
      return CGAL::compare(sqdistance(ap, h1), sqdistance(ap, h2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      // if (cv.is_segment()) {
      //   p = k.construct_midpoint_2_object()(cv.left(), cv.right());
      // }
      // else if (cv.is_ray()) {
      //   p = k.construct_point_on_2_object()(cv.ray(), 1);
      // }
      // else {
      //   CGAL_assertion(cv.is_line());
      //   p = k.construct_point_on_2_object()(cv.line(), 1);
      // }

      /* compare using the middle point */
      Alg_point_2 p = construct_middle_point(cv);
      return this->operator()(p, h1, h2);
    }

    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      /* if the two unbounded surfaces do not intersect, then they must
       * represent the same segment's distance function */
      return CGAL::EQUAL; // they are literally the same surface
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
      return CGAL::LARGER; //TODO fake
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
      return CGAL::LARGER; //TODO fake
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }

}; // class L2_segment_voronoi_traits_2
} // namespace CGAL

#endif // CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
