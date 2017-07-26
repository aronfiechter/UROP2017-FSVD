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

/* to compute bisector of two segments */
#include <CGAL/ch_akl_toussaint.h> // convex hull

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
  typedef typename CGAL::Polygon_2<Rat_kernel>        Rat_polygon_2;
  typedef typename Rat_polygon_2::Edge_const_iterator Edge_iterator;

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

  /* Construct a parabolic arc on the parabola with a directrix that is the line
   * that supports the segment seg and a focus that is the point f. The arc goes
   * from points p1 to p2.
   * Precondition: p1 and p2 are on the parabola */
  static Curve_2 construct_parabolic_arc(Rat_segment_2 seg, Rat_point_2 f,
    Point_2 p1, Point_2 p2) {

    /* get supporting_line of seg */
    Rat_line_2 directrix = seg.supporting_line();
    Rat_point_2 focus = f;

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

    /* construct the curve using the parameters and the endpoints; the
     * orientation is clockwise becuase probably it doesn't matter at all */
    Curve_2 arc(r, s, t, u, v, w, CGAL::CLOCKWISE, p1, p2);

    CGAL_assertion(arc.is_valid()); // valid arc
    CGAL_assertion(4 * r * s - CGAL::square(t) == 0); // curve is a parabola
    return arc;
  }

  /* Convert the Curve_2 cv into multiple X_monotone_curve_2 using the provided
   * make_x_monotone function. Store the results into the list x_mono_curves.
   * Precondition: cv is a valid curve. */
  static void make_curve_2_into_many_x_monotone_curve_2(Curve_2& cv,
    std::vector<X_monotone_curve_2>& x_mono_curves) {

    /* instantiate traits, we need the provided function */
    C_traits_2 c_traits;
    typename C_traits_2::Make_x_monotone_2 make_x_monotone =
      c_traits.make_x_monotone_2_object();

    /* call the provided function */
    std::vector<CGAL::Object> pre_x_mono_curves;
    make_x_monotone(cv, std::back_inserter(pre_x_mono_curves));

    /* cast all CGAL::Objects into X_monotone_segment_2 and add to list */
    for(size_t i = 0; i < pre_x_mono_curves.size(); i++ ) {
      X_monotone_curve_2 curr;
      bool check = CGAL::assign(curr, pre_x_mono_curves[i]);
      assert(check); CGAL_USE(check);
      x_mono_curves.push_back(curr);
    }
  }

  /* Check if a given segment called edge connects two the other segments s1 and
   * s2 by any of their endpoints.
   * Return true if this is the case, return false if edge is actyally just
   * connecting s1's endpoints (or s2's). We cannot just check for equality
   * because edge could be just the same as one segment but in the other
   * direction. */
  static bool edge_connects_segments(Rat_segment_2 edge, Rat_segment_2 s1,
    Rat_segment_2 s2) {

    /* create a copy of edge but in the other direction, then check equality for
     * both versions of edge */
    Rat_segment_2 rev_edge(edge.target(), edge.source());
    return !(edge == s1 || edge == s2 || rev_edge == s1 || rev_edge == s2);
  }



private:
  //TODO remove if unused
  // enum Rel_position {
  //   NO_INFLUENCE,
  //   PARTIAL_INFLUENCE,
  //   COMPLETE_INFLUENCE
  // };
  //
  // /* Determine the relative position of two segments in R_2. */
  // static Rel_position relative_position(Rat_segment_2 s1, Rat_segment_2 s2) {
  //   //TODO fake
  //   return COMPLETE_INFLUENCE;
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
      Rat_kernel rat_kernel;
      /* if the two segments are the same, their distance function is the same,
       * so there is no intersection */
      if (s1 == s2) {
        return o;
      }
      /* if the two segments do not intersect, construct the bisector starting
       * from one unbounded
       */
      else if (!CGAL::do_intersect(s1, s2)) {
        /* first of all compute the two unbounded edges of the bisector, which
         * have to be saved as two very long segments. To do this, first compute
         * the convex hull of the endpoints of the segments. The two pairs of
         * vertices of the hull that are not of the same segment are the pairs
         * of which the bisector lines contain the two unbouded rays that are
         * the unbounded rays of the plane bisector of the two segments */

        /* compute hull of endpoints */
        std::list<Rat_point_2> ch_points;
        std::list<Rat_point_2> points = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };
        CGAL::ch_akl_toussaint(
          points.begin(), points.end(),
          std::back_insert_iterator<std::list<Rat_point_2>>(ch_points)
        );

        /* make a polygon out of the hull points, iterate over vetrices to find
         * pairs to make rays, directed towards outside of polygon */
        Rat_polygon_2 ch_polygon(ch_points.begin(), ch_points.end());
        CGAL_assertion(ch_polygon.is_convex()); // it is a hull
        for ( // for all edges
          Edge_iterator eit = ch_polygon.edges_begin();
          eit != ch_polygon.edges_end();
          ++eit
        ) {
          if (edge_connects_segments(*eit, s1, s2)) {
            Rat_segment_2 seg = *eit;
            X_monotone_curve_2 curve_seg(seg);
            *o++ = CGAL::make_object(
              Intersection_curve(curve_seg, 0)
            ); //TODO this just draws on Ipe all connecting edges of the hull
          }
        }

        // Curve_2 par_arc = construct_parabolic_arc(s1, s2.source(), i1, i2);

        /* critical, if the curve is not valid, abort immediately */
        // if (!c2.is_valid()) return o;
        // *o++ = CGAL::make_object(
        //   Intersection_curve(c2, 0)
        // );

        // *o++ = CGAL::make_object(
        //   Intersection_curve(X_monotone_curve_2(Rat_segment_2(s1.source(), s2.source())), 0)
        // );

        return o;
      }
      /* if instead they do intersect, assert it, then proceed to computing then
       * intersection in this case */
      else {
        CGAL_assertion(CGAL::do_intersect(s1, s2)); // they HAVE to intersect
        return o; //TODO fake
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
