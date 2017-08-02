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
  typedef Kernel_                                     Kernel;
  typedef Conic_traits_2                              C_traits_2;
  typedef L2_segment_voronoi_traits_2<C_traits_2, Kernel> Self;

  typedef typename C_traits_2::Point_2                Point_2;
  typedef typename C_traits_2::Curve_2                Curve_2;
  typedef typename C_traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename C_traits_2::Multiplicity           Multiplicity;

  typedef typename C_traits_2::Rat_kernel             Rat_kernel;
  typedef typename C_traits_2::Alg_kernel             Alg_kernel;
  typedef typename C_traits_2::Nt_traits              Nt_traits;

  typedef typename Rat_kernel::FT                     Rational;
  typedef typename Rat_kernel::Point_2                Rat_point_2;
  typedef typename Rat_kernel::Segment_2              Rat_segment_2;
  typedef typename Rat_kernel::Line_2                 Rat_line_2;
  typedef typename Rat_kernel::Ray_2                  Rat_ray_2;
  typedef typename Rat_kernel::Vector_2               Rat_vector_2;
  typedef typename Rat_kernel::Direction_2            Rat_direction_2;
  typedef typename CGAL::Polygon_2<Rat_kernel>        Rat_polygon_2;
  typedef typename Rat_polygon_2::Edge_const_iterator Edge_iterator;

  typedef typename Alg_kernel::FT                     Algebraic;
  typedef typename Alg_kernel::Point_2                Alg_point_2;
  typedef typename Alg_kernel::Segment_2              Alg_segment_2;
  typedef typename Alg_kernel::Line_2                 Alg_line_2;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Rat_segment_2                               Surface_3;
  typedef Surface_3                                   Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity> Intersection_curve;

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

  /* Given a bisector finds the point that is the furthest intersection
   * (following the direction of the bisector) of the bisector with the four
   * lines saved in delimiters */
  static Rat_point_2 find_unbounded_ray_start_point(
    Rat_line_2 bisector,
    std::pair<
      std::pair<Rat_line_2, Rat_line_2>,
      std::pair<Rat_line_2, Rat_line_2>
    > delimiters
  ) {

    /* get the four intersection points, add them to two lists, sort one by x
     * and the other one by y */
    Rat_point_2 p1; // intersection between bisector and delimiters[0][0]
    Rat_point_2 p2; // intersection between bisector and delimiters[0][1]
    Rat_point_2 p3; // intersection between bisector and delimiters[1][0]
    Rat_point_2 p4; // intersection between bisector and delimiters[1][1]
    CGAL::assign(p1, CGAL::intersection(bisector, delimiters.first.first));
    CGAL::assign(p2, CGAL::intersection(bisector, delimiters.first.second));
    CGAL::assign(p3, CGAL::intersection(bisector, delimiters.second.first));
    CGAL::assign(p4, CGAL::intersection(bisector, delimiters.second.second));
    std::vector<Rat_point_2> intersection_x = { p1, p2, p3, p4 };
    std::vector<Rat_point_2> intersection_y = { p1, p2, p3, p4 };
    std::sort(intersection_x.begin(), intersection_x.end(),
      [](Rat_point_2 a, Rat_point_2 b) {
      return a.x() < b.x();
    });
    std::sort(intersection_y.begin(), intersection_y.end(),
      [](Rat_point_2 a, Rat_point_2 b) {
      return a.y() < b.y();
    });

    /* find the farthest point according to the direction of bisector */
    Rat_direction_2 dir = bisector.direction();
    if (dir.dx() == 0) {
      return (dir.dy() > 0) ? intersection_y.back() : intersection_x.front();
    }
    else {
      return (dir.dx() > 0) ? intersection_x.back() : intersection_x.front();
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
      /* the surfaces we are considering are distance functions of line
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
  private:
      Rational UNBOUNDED_RAY_LENGTH = Rational(10, 1);

  public:
    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {
      Rat_kernel rat_kernel;
      /* if the two segments are the same (also if one is just the other but
       * reversed), their distance function is the same, so there is no
       * intersection */
      if (s1 == s2 || s1 == Rat_segment_2(s2.target(), s2.source())) {
        return o;
      }
      /* if the two segments are not the same, compute all parts of their plane
       * bisector */
      else {
        /* first of all, for each segment create the two lines that divide the
         * plane in three areas: one of all points closest to the inner part of
         * the segment, the other two of all points closest to the two endpoints
         * of the segment.
         * The lines are saved with an orientation such that they both have the
         * inner part of the segment on their right side (negative side).
         * Note: a vector constructed using a segment is oriented from source to
         * target of that segment, so to build a line such that the segment lies
         * on the right side of it, we need to use:
         * - the source of the segment and as direction the vector oriented 90
         *   degrees counterclockwise from the segment vector
         * - the target of the segment but as direction the vector oriented 90
         *   degrees clockwise. */
        std::pair<
          std::pair<Rat_line_2, Rat_line_2>,
          std::pair<Rat_line_2, Rat_line_2>
        > delimiter_lines = {
          {
            Rat_line_2(
              s1.source(),
              Rat_vector_2(s1).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s1.target(),
              Rat_vector_2(s1).perpendicular(CGAL::CLOCKWISE)
            )
          },
          {
            Rat_line_2(
              s2.source(),
              Rat_vector_2(s2).perpendicular(CGAL::COUNTERCLOCKWISE)
            ),
            Rat_line_2(
              s2.target(),
              Rat_vector_2(s2).perpendicular(CGAL::CLOCKWISE)
            )
          }
        };
        /* also save the lines in a vector for convenience */
        std::vector<Rat_line_2> delimiter_lines_vector = {
          delimiter_lines.first.first, delimiter_lines.first.second,
          delimiter_lines.second.first, delimiter_lines.second.second
        };

        /* then compute the 2 or 4 unbounded edges of the bisector.
         * To do this, first compute the convex hull of the endpoints of the
         * segments. The pairs of vertices of the hull that are not of the same
         * segment are the pairs of which the bisector lines contain the
         * unbounded rays that are the unbounded rays of the plane bisector of
         * the two segments */

        /* compute hull of endpoints; the hull points will be stored inside
         * ch_points in counterclockwise order */
        std::list<Rat_point_2> ch_points;
        std::list<Rat_point_2> points = {
          s1.source(), s1.target(), s2.source(), s2.target()
        };
        CGAL::ch_akl_toussaint(
          points.begin(), points.end(),
          std::back_insert_iterator<std::list<Rat_point_2>>(ch_points)
        );

        /* make a polygon out of the hull points, iterate over vertices to find
         * pairs to make rays, directed towards outside of polygon */
        Rat_polygon_2 ch_polygon(ch_points.begin(), ch_points.end());
        CGAL_assertion(ch_polygon.is_convex()); // it is a hull
        CGAL_assertion(ch_polygon.area() >= 0); // it is counterclockwise

        /* list to save starting points of unbounded rays */
        std::list<Rat_point_2> ray_start_points;

        for ( // for all edges
          Edge_iterator eit = ch_polygon.edges_begin();
          eit != ch_polygon.edges_end();
          ++eit
        ) {
          if (edge_connects_segments(*eit, s1, s2)) { // if it's not s1 or s2
            /* create line that bisects the segment, orient it outside */
            Rat_line_2 bisector_line = CGAL::bisector(
              eit->target(), eit->source()
            );
            Rat_point_2 start_point = find_unbounded_ray_start_point(
              bisector_line, delimiter_lines
            );
            ray_start_points.push_back(start_point);
            Rat_point_2 end_point = start_point
              + UNBOUNDED_RAY_LENGTH * bisector_line.direction().vector();

            /* make very long segment to represent an unbounded ray, so that it
             * can be saved as an X_monotone_curve_2, because the Conic_traits
             * require that curves are bounded */
            Rat_segment_2 seg(start_point, end_point);
            X_monotone_curve_2 curve_seg(seg);
            *o++ = CGAL::make_object(
              Intersection_curve(curve_seg, 0)
            );
          }
        }

        /* if the two segments do NOT intersect, construct the bisector starting
         * from one unbounded edge, finding the correct intersection points
         * using the delimiter_lines.
         * In this case, the ray start points should be only two. */
        if (!CGAL::do_intersect(s1, s2)) { // segments do not intersect
          CGAL_assertion(ray_start_points.size() == 2);

          /* start from one point, find intersection line */
          Rat_point_2 start = ray_start_points.front();
          Rat_line_2 intersecting_delimiter;
          bool found;
          std::for_each(
            delimiter_lines_vector.begin(),
            delimiter_lines_vector.end(),
            [&intersecting_delimiter, &start, &found] (Rat_line_2 delimiter) {
            if (delimiter.has_on(start)) {
              found = true;
              intersecting_delimiter = delimiter;
            }
          });
          CGAL_assertion(found);
          // TODO revise this horrendous thing, could get lines already
          // by returning a pair from find_unbounded_ray_start_point

          

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
        } // end of segments do not intersect
        /* if instead they do intersect, assert it, then proceed to computing
         * the bisector in this case.
         * In this case, the ray start points should be four. */
        else {
          CGAL_assertion(CGAL::do_intersect(s1, s2)); // they HAVE to intersect
          CGAL_assertion(ray_start_points.size() == 4);
          return o;
        } // end of segments intersect

      } // end of segments are not the same
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
