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
  typedef Rat_segment_2                             Surface_3;
  typedef Surface_3                                 Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity>    Intersection_curve;

  /* Returns the squared distance between two points in L2 metric. */
  static Algebraic sqdistance(const Point_2& p1, const Point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static Algebraic sqdistance(const Point_2& p, const Rat_segment_2& s) {
    /* find projection of p on supporting line of s */
    Alg_segment_2 alg_s(
      Alg_point_2(s.source().x(), s.source().y()),
      Alg_point_2(s.target().x(), s.target().y())
    );
    Alg_line_2 l = alg_s.supporting_line();
    Alg_point_2 proj = l.projection(p);

    /* if the projection is on s, the distance is d(p,proj) */
    if (alg_s.has_on(proj)) {
      return sqdistance(p, proj);
    }
    /* otherwise, the distance is min(d(p,s1),d(p,s2)), where s1 and s2 are the
    endpoints of the segment s */
    else {
      return CGAL::min(
        sqdistance(p, Alg_point_2(s.source().x(), s.source().y())),
        sqdistance(p, Alg_point_2(s.target().x(), s.target().y()))
      );
    }
  }

  /* Construct a point in the middle of the curve cv. This function is copied
   * from Env_sphere_traits_3.h */
  static Point_2 construct_middle_point(const X_monotone_curve_2& cv) {
    /* get the x-value of the middle point */
    Alg_kernel k;
    Alg_point_2 mid_x = k.construct_midpoint_2_object()(
      cv.source(),
      cv.target()
    );

    /* if cv is vertical, it is just a segment */
    if (cv.is_vertical()) return Point_2(mid_x);
    /* otherwise take the point with the same x coordinate but on cv */
    else return Point_2(cv.point_at_x(mid_x));
  }

  /* Construct a parabolic arc that is the bisector between one endpoint of one
   * segment and the inner part of the other. //TODO fake */
  static Curve_2 construct_parabolic_arc(Rat_segment_2 s1, Rat_segment_2 s2) {
    Rat_line_2 directrix = s1.supporting_line();
    Rat_point_2 focus = s2.target();

    Rational a = directrix.a();
    Rational b = directrix.b();
    Rational c = directrix.c();
    Rational f_x = focus.x();
    Rational f_y = focus.y();
    Rational NEG2 = Rational(-2);

    Rational r = CGAL::square(b);
    Rational s = CGAL::square(a);
    Rational t = NEG2 * a * b;
    Rational u =
      NEG2 * a * c +
      NEG2 * CGAL::square(a) * f_x +
      NEG2 * CGAL::square(b) * f_x
    ;
    Rational v =
      NEG2 * b * c +
      NEG2 * CGAL::square(a) * f_y +
      NEG2 * CGAL::square(b) * f_y
    ;
    Rational w =
      CGAL::square(a) * CGAL::square(f_x) +
      CGAL::square(a) * CGAL::square(f_y) +
      CGAL::square(b) * CGAL::square(f_x) +
      CGAL::square(b) * CGAL::square(f_y) -
      CGAL::square(c)
    ;

    Curve_2 arc(r, s, t, u, v, w);

    CGAL_assertion(arc.is_valid()); // valid arc
    CGAL_assertion(4 * r * s - CGAL::square(t) == 0); // curve is a parabola
    return arc;
  }



private:
  enum Rel_position {
    NO_INFLUENCE,
    PARTIAL_INFLUENCE,
    COMPLETE_INFLUENCE
  };

  /* Determine the relative position of two segments in R_2. */
  static Rel_position relative_position(Rat_segment_2 s1, Rat_segment_2 s2) {
    //TODO fake
    return COMPLETE_INFLUENCE;
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



  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {
      Rat_kernel rat_kernel;
      /* if the two segments are the same, their distance function is the same,
       * so there is no intersection */
      if (s1 == s2) {
        return o;
      }
      /* if they do not intersect, there are three main cases:
       * - the areas of the points closest to the inner parts of both segments
       *   do not intersect at all (the bisector has 5 parts)
       * - the area of the points closest to the inner part of one segment
       *   intersects the other segment (the bisector has 5 parts)
       * - the area of the points closest to the inner part of one segment
       *   completely contains the other segment (the bisector has 7 parts)
       *
       * each one of the cases has many subcases, depending also on the distance
       * between the segments and their length
       */
      else {
        switch (relative_position(s1, s2)) {
          case NO_INFLUENCE: {
            //TODO not implemented
            break;
          }

          case PARTIAL_INFLUENCE: {
            //TODO not implemented
            break;
          }

          case COMPLETE_INFLUENCE: {
            //TODO not implemented
            construct_parabolic_arc(s1, s2);
            break;
          }

          default: {
            break;
          }
        }

        // Curve_2 c1(
        //   1, 1, 1, 4, -1, 0, CGAL::COUNTERCLOCKWISE,
        //   Point_2(s1.source().x(), s1.source().y()),
        //   Point_2(s2.target().x(), s2.target().y())
        // );
        // CGAL_assertion(c1.is_valid());
        // X_monotone_curve_2 mc1(c1);
        // *o++ = CGAL::make_object(
        //   Intersection_curve(mc1, 0)
        // );
        C_traits_2 c_traits;
        typename C_traits_2::Make_x_monotone_2 make_x_monotone = c_traits.make_x_monotone_2_object();

        Curve_2 c2(
          Rat_point_2(s1.source()),
          Rat_point_2(s2.source()),
          Rat_point_2(s1.target())
        );

        /* critical, if the curve is not valid, abort immediately */
        if (!c2.is_valid()) return o;

        // std::vector<X_monotone_curve_2> segs;
        std::vector<CGAL::Object> pre_segs;
        make_x_monotone(c2, std::back_inserter(pre_segs));

        for (size_t i = 0; i < pre_segs.size(); i++ ) {
          X_monotone_curve_2 curr;
          bool check = CGAL::assign(curr, pre_segs[i]);
          assert(check); CGAL_USE(check);
          // segs.push_back(curr);

          *o++ = CGAL::make_object(
            Intersection_curve(curr, 0)
          );
        }

        // *o++ = CGAL::make_object(
        //   Intersection_curve(mc2, 0)
        // );
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
      return CGAL::compare(sqdistance(p, h1), sqdistance(p, h2));
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
      Point_2 p = construct_middle_point(cv);
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
