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

/* to convert from Alg to Rat and viceversa */
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Cartesian.h>

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

  typedef Rational                                    RT;
  typedef Algebraic                                   AT;

  /* Converters */
  typedef CGAL::Cartesian_converter<Alg_kernel, Rat_kernel> AK_to_RK;
  typedef CGAL::Cartesian_converter<Rat_kernel, Alg_kernel> RK_to_AK;
  typedef CGAL::Cartesian<double>                     D_kernel;
  typedef typename D_kernel::Point_2                  D_point_2;
  typedef CGAL::Cartesian_converter<Alg_kernel, D_kernel> AK_to_DK;
  typedef CGAL::Cartesian_converter<D_kernel, Rat_kernel> DK_to_RK;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Rat_segment_2                               Surface_3;
  typedef Surface_3                                   Xy_monotone_surface_3;

protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity> Intersection_curve;

private:
  enum Seg_endpoint {
    S1_SOURCE,
    S1_TARGET,
    S2_SOURCE,
    S2_TARGET,
  };
  typedef typename std::pair<Rat_ray_2, Seg_endpoint> Ray_info;

  typedef typename std::pair<
    std::pair<Rat_line_2, Rat_line_2>,
    std::pair<Rat_line_2, Rat_line_2>
  >                                                   Delimiter_lines;

  enum Bisector_type {
    PARABOLIC_ARC,
    SUPP_LINE_BISECTOR,
    ENDPOINT_BISECTOR
  };

  class Parabola {

  private:

    /* Fields */

    /* coefficients of equation: rx^2 + sy^2 + txy + ux + vy + w = 0 */
    RT _r; RT _s; RT _t; RT _u; RT _v; RT _w;

    /* generators of parabola, direction of parabola is the same of directrix */
    Rat_line_2 _directrix;
    Rat_point_2 _focus;

    /* Private constructor to save equation coefficients */
    Parabola(RT __r, RT __s, RT __t, RT __u, RT __v, RT __w)
      : _r(__r), _s(__s), _t(__t), _u(__u), _v(__v), _w(__w) {
      CGAL_assertion(CGAL::square(_t) - 4 * _r * _s == 0); // curve is parabola
    }

  public:

    /* Empty constructor */
    Parabola() {}

    /* Construct using directrix and focus. Details on the computation of the
     * forumla can be found in doc/parabola.pdf */
    Parabola(Rat_line_2 directrix, Rat_point_2 focus)
      : _directrix(directrix), _focus(focus) {

      RT a = directrix.a();
      RT b = directrix.b();
      RT c = directrix.c();
      RT f_x = focus.x();
      RT f_y = focus.y();
      RT TWO = RT(2);

      RT r = -CGAL::square(b);
      RT s = -CGAL::square(a);
      RT t = TWO * a * b;
      RT u =
        TWO * a * c +
        TWO * CGAL::square(a) * f_x +
        TWO * CGAL::square(b) * f_x
      ;
      RT v =
        TWO * b * c +
        TWO * CGAL::square(a) * f_y +
        TWO * CGAL::square(b) * f_y
      ;
      RT w =
        CGAL::square(c) -
        CGAL::square(a) * CGAL::square(f_x) -
        CGAL::square(a) * CGAL::square(f_y) -
        CGAL::square(b) * CGAL::square(f_x) -
        CGAL::square(b) * CGAL::square(f_y)
      ;

      /* finish constructing using the private coefficients constructor */
      Parabola(r, s, t, u, v, w);
    }

    /* Getters */
    RT r() { return _r; }
    RT s() { return _s; }
    RT t() { return _t; }
    RT u() { return _u; }
    RT v() { return _v; }
    RT w() { return _w; }
    Rat_line_2 directrix() { return _directrix; }
    Rat_point_2 focus() { return _focus; }

    /* Methods */

    /* Evaluate the equation of the parabola rx^2 + sy^2 + txy + ux + vy + w = 0
     * using the x and y of the point. If the result is 0, the point is on the
     * parabola, if it is positive the point lies on the positive side, if it
     * is negative the point lies on the negative side. */
    Algebraic evaluate(Point_2 point) {
      Algebraic x = point.x();
      Algebraic y = point.y();
      return Algebraic(
        this->r() * CGAL::square(x) +
        this->s() * CGAL::square(x) +
        this->t() * x * y +
        this->u() * x +
        this->v() * y +
        this->w()
      );
    }

    /* Check if a given point lies on the parabola by checking if the values of
     * x and y (the point's coordinates) satisfy the equation of the parabola */
    bool has_on(Point_2 point) {
      return CGAL::is_zero(this->evaluate(point));
    }
    /* Check if a given point lies on the positive side of the parabola. The
     * positive side is the one on the left when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_positive_side(Point_2 point) {
      return CGAL::is_positive(this->evaluate(point));
    }
    /* Check if a given point lies on the negative side of the parabola. The
     * negative side is the one on the right when traveling on the curve in the
     * same direction as the directrix (by construction they have the "same"
     * oriented side) */
    bool has_on_negative_side(Point_2 point) {
      return CGAL::is_negative(this->evaluate(point));
    }

    /* Save into the OutputIterator o the intersection(s) of the parabola with
     * a given line l. The type of o must be Alg_point_2.
     * Return a past the end iterator o. */
    template <class OutputIterator>
    OutputIterator get_intersections(Rat_line_2 line, OutputIterator o) {
      /* equation of line:      ax + by + c = 0
       * equation of parabola:  rx^2 + sy^2 + txy + ux + vy + w = 0
       * we can find intersections by substituting line in parabola; described
       * in detail in docs/parabola.pdf, verified using wxMaxima */
      RT a = line.a();
      RT b = line.b();
      RT c = line.c();

      /* convert line to algebraic, get nt_traits to solve quadratic equation */
      RK_to_AK to_alg;
      Alg_line_2 alg_line = to_alg(line);
      Nt_traits nt_traits;

      /* in this case the intersection is simpler, since we can substitute the x
       * in the parabola equation with just:    x = -c/a
       * We get a quadratic equation in y, which we can solve using CGAL */
      if (b == 0) {
        /* the quadratic equation in y is:
         *                                    s y^2
         *                     + (v - (ct / a)) y
         *      + ((rc^2 / a^2) - (cu / a) + w)
         *                                          = 0
         */
        RT EQ_A = this->s();
        RT EQ_B = this->v() - ((c * this->t()) / a);
        RT EQ_C = ((this->r() * CGAL::square(c)) / CGAL::square(a)) -
          ((c * this->u()) / a) +
          this->w()
        ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting y, find the corresponding x, and add a point to
         * the OutputIterator o */
        Algebraic  ys[2];
        Algebraic * ys_end;
        int n_ys;
        ys_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, ys);
        n_ys = ys_end - ys;

        /* if no intersections return */
        if (n_ys == 0) {
          return o;
        }
        /* else find xs for all ys, add points to iterator */
        else while (--n_ys >= 0) {
          Algebraic current_y = ys[n_ys];
          Algebraic corresponding_x = alg_line.x_at_y(current_y);
          *o++ = Alg_point_2(corresponding_x, current_y);
        }

        return o; // already one past the end, post-incremented when adding
      }
      /* in the general case we substitute the y in the parabola equation with
       * the value:                             y = -c/b + -ax/b
       * We get a quadratic equation in x, which we can solve using CGAL */
      else {
        /* the quadratic equation in x is:
         *                (r + (sa^2 / b^2) - (at / b)) x^2
         *   + ((2acs / b^2) - (ct / b) - (av / b) + u) x
         *              + ((sc^2 / b^2) - (cv / b) + w)
         *                                                  = 0
         */
        RT EQ_A = this->r() +
          (this->s() * CGAL::square(a) / CGAL::square(b)) -
          (a * this->t() / b)
        ;
        RT EQ_B = (2 * a * c * this->s() / CGAL::square(b)) -
          (c * this->t() / b) -
          (a * this->v() / b) +
          this->u()
        ;
        RT EQ_C = (this->s() * CGAL::square(c) / CGAL::square(b)) -
          (c * this->v() / b) +
          this->w()
        ;

        /* to store the 0, 1, or 2 results, indicating the intersections.
         * For all resulting x, find the corresponding y, and add a point to
         * the OutputIterator o */
        Algebraic  xs[2];
        Algebraic * xs_end;
        int n_xs;
        xs_end = nt_traits.solve_quadratic_equation(EQ_A, EQ_B, EQ_C, xs);
        n_xs = xs_end - xs;

        /* if no intersections return */
        if (n_xs == 0) {
          return o;
        }
        /* else find xs for all xs, add points to iterator */
        else while (--n_xs >= 0) {
          Algebraic current_x = xs[n_xs];
          Algebraic respective_y = alg_line.y_at_x(current_x);
          *o++ = Alg_point_2(respective_y, current_x);
        }

        return o; // already one past the end, post-incremented when adding
      }
    }

    /* Construct a parabolic arc on the parabola from point p1 to point p2.
     * Precondition (checked): p1 and p2 are on the parabola */
    Curve_2 construct_parabolic_arc(Point_2 p1, Point_2 p2) {
      /* check precondition: both points lie on the parabola */
      CGAL_assertion(this->has_on(p1));
      CGAL_assertion(this->has_on(p2));

      /* construct the curve using the parameters and the endpoints */
      Curve_2 arc(_r,_s,_t,_u,_v,_w, CGAL::CLOCKWISE, p1, p2); //TODO ORIENTATION

      CGAL_assertion(arc.is_valid()); // valid arc
      return arc;
    }
  };

  /* Returns the squared distance between two points in L2 metric. */
  static Algebraic sqdistance(const Point_2& p1, const Point_2& p2) {
    return CGAL::squared_distance(p1, p2);
  }
  static Rational sqdistance(const Rat_point_2& p1, const Rat_point_2& p2) {
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

  /* Given a point p and a list of points, return the closest point */
  static Rat_point_2 closest_pt(Rat_point_2 p, std::list<Rat_point_2> points) {
    Rat_point_2 result;
    RT smaller_sqdistance = -1;
    for (auto& q : points) {
      RT sqdist_pq = sqdistance(p, q);
      if (smaller_sqdistance < 0 || smaller_sqdistance > sqdist_pq) {
        result = q;
        smaller_sqdistance = sqdist_pq;
      }
    }

    CGAL_assertion(smaller_sqdistance >= 0);
    return result;
  }

  /* Construct a point in the middle of the curve cv. This function is copied
   * from Env_sphere_traits_3.h */
  static Point_2 construct_middle_point(const X_monotone_curve_2& cv) {
    /* get the x-value of the middle point */
    Alg_point_2 mid_x = CGAL::midpoint(cv.source(),cv.target());

    /* if cv is vertical, it is just a segment */
    if (cv.is_vertical()) return Point_2(mid_x);
    /* otherwise take the point with the same x coordinate but on cv */
    else return Point_2(cv.point_at_x(mid_x));
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
   * lines saved in delimiters.
   * Return a pair with the unbounded ray and the segment endpoint whose
   * orthogonal delimiter intersects the ray's source. */
  static Ray_info find_unbounded_ray(
    Rat_line_2 bisector,
    Delimiter_lines delimiters
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
    std::vector<Ray_info> intersections_x = {
      std::make_pair(Rat_ray_2(p1, bisector.direction()), S1_SOURCE),
      std::make_pair(Rat_ray_2(p2, bisector.direction()), S1_TARGET),
      std::make_pair(Rat_ray_2(p3, bisector.direction()), S2_SOURCE),
      std::make_pair(Rat_ray_2(p4, bisector.direction()), S1_TARGET)
    };
    std::vector<Ray_info> intersections_y;
    std::copy(intersections_x.begin(), intersections_x.end(),
              std::back_inserter(intersections_y));
    std::sort(intersections_x.begin(), intersections_x.end(),
      [](Ray_info a, Ray_info b) {
      return a.first.source().x() < b.first.source().x();
    });
    std::sort(intersections_y.begin(), intersections_y.end(),
      [](Ray_info a, Ray_info b) {
      return a.first.source().y() < b.first.source().y();
    });

    /* find the farthest point according to the direction of bisector */
    Rat_direction_2 dir = bisector.direction();
    if (dir.dx() == 0) {
      return (dir.dy() > 0) ? intersections_y.back() : intersections_x.front();
    }
    else {
      return (dir.dx() > 0) ? intersections_x.back() : intersections_x.front();
    }
  }


  /* Given a direction and a start point finds the point that is the first
   * intersection (after this start point with the four lines saved in the
   * vector of delimiter lines.
   * Return the found intersection point. */
  static Rat_point_2 find_next_intersection(
    Rat_direction_2 direction,
    Rat_point_2 start_pt,
    std::vector<Rat_line_2> delimiters
  ) {
    /* list to store intersections */
    std::list<Rat_point_2> intersections;

    /* for each delimiter add the intersection with ray, if it's not start_pt */
    Rat_ray_2 ray(start_pt, direction);
    for (auto& delimiter : delimiters) {
      if (CGAL::do_intersect(delimiter, ray)) {
        Rat_point_2 intersection;
        CGAL::assign(intersection, CGAL::intersection(delimiter, ray));
        if (intersection != start_pt) intersections.push_back(intersection);
      }
    }

    /* all intersections are in the correct direction because we used a ray
     * starting from start_pt, so return the closest one */
    return closest_pt(start_pt, intersections);
  }

  /* Determine the position of the point p relative to the segments s1 and s2.
   * We do not have the segments though: we have four lines, each orthogonal to
   * one endpoint of one of the two segments. Every line is oriented so to have
   * the inner part of the segment on their negative side (right side).
   * There are three main cases (described by enum Bisector_type):
   * - PARABOLIC_ARC: when p is closer to one segment's inner part and to one of
   *   the other segment's endpoints. In this case, save the supporting_line of
   *   the first segment and the endpoint of the second segment.
   * - SUPP_LINE_BISECTOR: when p is closer to both inner parts of both the two
   *   segments. In this case save the two supporting_lines of the two segments.
   * - ENDPOINT_BISECTOR: when p is closer to two endpoints of the two segments.
   *   In this case save those two endpoints.
   * In all three cases we save in o1 the correct endpoint or supporting_line of
   * s1, and in o2 the same for s2.
   */
  static Bisector_type find_position(
    Alg_point_2 p,
    Delimiter_lines delimiter_lines,
    Rat_segment_2 s1,
    Rat_segment_2 s2,
    Object& o1,  // to store [directrix1 or focus1]/line1/point1
    Object& o2   // to store [directrix2 or focus2]/line2/point2
  ) {
    /* convert point p from alg to rational */
    AK_to_DK to_dbl;
    DK_to_RK to_rat;
    Rat_point_2 rp = to_rat(to_dbl(p)); //TODO revise to find better solution

    /* assume point is not on any delimiter, consider all other cases. To do so,
     * first determine what must be stored in o1, then in o2. Save in two flags
     * information about the case. In the end, determine the case. */
    bool o1_is_line = false;
    bool o2_is_line = false;

    /* determine o1 */
    if (delimiter_lines.first.first.has_on_positive_side(rp)) {
      o1 = CGAL::make_object(s1.source());
    }
    else if (delimiter_lines.first.second.has_on_positive_side(rp)) {
      o1 = CGAL::make_object(s1.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.first.first.has_on_negative_side(rp)
        &&
        delimiter_lines.first.second.has_on_negative_side(rp)
      );
      o1 = CGAL::make_object(s1.supporting_line());
      o1_is_line = true;
    }

    /* determine o2 */
    if (delimiter_lines.second.first.has_on_positive_side(rp)) {
      o2 = CGAL::make_object(s2.source());
    }
    else if (delimiter_lines.second.second.has_on_positive_side(rp)) {
      o2 = CGAL::make_object(s2.target());
    }
    else {
      CGAL_assertion(
        delimiter_lines.second.first.has_on_negative_side(rp)
        &&
        delimiter_lines.second.second.has_on_negative_side(rp)
      );
      o2 = CGAL::make_object(s2.supporting_line());
      o2_is_line = true;
    }

    /* determine case using flags */
    if (o1_is_line) {
      return (o2_is_line) ? SUPP_LINE_BISECTOR : PARABOLIC_ARC;
    }
    else {
      return (o2_is_line) ? PARABOLIC_ARC : ENDPOINT_BISECTOR;
    }
  }

  /* Given a line and a direction determine whether the line is oriented in that
   * genreal direction, that is in a range of [-90˚, 90˚] around the Given
   * direction. */
  static bool generally_same_direction(Rat_line_2 line, Rat_direction_2 dir) {
    return line.direction().counterclockwise_in_between(
      dir.vector().perpendicular(CGAL::CLOCKWISE).direction(),
      dir.vector().perpendicular(CGAL::COUNTERCLOCKWISE).direction()
    ); //TODO what if they are perpendicular? which case is it?
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
        Delimiter_lines delimiter_lines = {
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
        /* also save segment endpoints "generating" these lines */
        std::vector<Rat_point_2> segment_endpoints = {
          s1.source(), s1.target(), s2.source(), s2.target()
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

        /* list to save starting points of unbounded rays, together with an
         * indication of which segment endpoint generates the orthogonal line
         * that caused the intersection */
        std::list<Ray_info> ray_info_list;

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
            /* find "farthest" intersection with delimiters, ray starts there */
            Ray_info ray_info = find_unbounded_ray(
              bisector_line, delimiter_lines
            );
            ray_info_list.push_back(ray_info);

            /* make very long segment to represent an unbounded ray, so that it
             * can be saved as an X_monotone_curve_2, because the Conic_traits
             * require that curves are bounded */
            Rat_point_2 start_point = ray_info.first.source();
            Rat_point_2 end_point = start_point
              + UNBOUNDED_RAY_LENGTH * ray_info.first.direction().vector();
            Rat_segment_2 seg(start_point, end_point);
            X_monotone_curve_2 curve_seg(seg);
            *o++ = CGAL::make_object(
              Intersection_curve(curve_seg, 0)
            );
          }
        }

        /* create converter functor to convert from Rational to Algebraic */
        RK_to_AK to_alg;

        /* if the two segments do NOT intersect, construct the bisector starting
         * from one unbounded edge, finding the correct intersection points
         * using the delimiter_lines.
         * In this case, the ray start points should be only two. */
        if (!CGAL::do_intersect(s1, s2)) { // segments do not intersect
          CGAL_assertion(ray_info_list.size() == 2);

          /* starting from the source of one unbounded ray and finishing at the
           * source of the other, compute the rest of the bisector, consisting
           * of:
           * - parabolic arcs: when we are in the "area of influence" of the
           *   interior of a segment and of one endpoint of the other
           * - segments: when we are in the "area of influence" of the interiors
           *   of the two segments or of two endpoints of the two segments */
          Ray_info start_ray_info = ray_info_list.front();
          Ray_info end_ray_info = ray_info_list.back();
          Rat_point_2 start_pt = start_ray_info.first.source();
          Rat_point_2 end_pt = end_ray_info.first.source();
          Rat_point_2 curr_pt = start_pt;
          Rat_direction_2 curr_direction = -(start_ray_info.first.direction());

          /* "walk" through the bisector to find all parts until every piece has
           * been created and added to the OutputIterator o */
          while (curr_pt != end_pt) {
            /* find next intersection with delimiter_lines when going in the
             * direction saved in "curr_direction", then find a middle point
             * between curr_pt and that intersection */
            Rat_point_2 approximate_next_intersection = find_next_intersection(
              curr_direction, curr_pt, delimiter_lines_vector
            );
            Rat_point_2 midpoint = CGAL::midpoint(
              curr_pt,
              approximate_next_intersection
            );

            /* determine where this middle point is relative to the two segments
             * s1 and s2, and create the correct piece of the bisector. The
             * objects o1 and o2 that are passed will store in the cases:
             * - PARABOLIC_ARC:       o1 = focus/directrix  o2 = focus/directrix
             * - SUPP_LINE_BISECTOR:  o1 = supp_line1,      o2 = supp_line2
             * - ENDPOINT_BISECTOR:   o1 = endpoint_1,      o2 = endpoint_2   */
            Object o1, o2;
            Curve_2 piece_of_bisector;
            Alg_point_2 alg_m_pt = to_alg(midpoint);
            switch (find_position(alg_m_pt, delimiter_lines, s1, s2, o1, o2)) {

              case PARABOLIC_ARC: {
                /* extract directrix and focus */
                Rat_line_2 directrix; Rat_point_2 focus;
                if (CGAL::assign(directrix, o1)) {
                  CGAL_assertion(CGAL::assign(focus, o2));
                }
                else {
                  CGAL_assertion(CGAL::assign(focus, o1));
                  CGAL_assertion(CGAL::assign(directrix, o2));
                }

                /* keep or invert directrix based on curr_direction */
                if (!generally_same_direction(directrix, curr_direction)) {
                  directrix = directrix.opposite();
                }

                /* create parabola */
                Parabola supporting_conic(directrix, focus);
                CGAL_assertion(supporting_conic.has_on(to_alg(curr_pt)));

                /* find actual next intersection of parabola */

                /* get parabolic arc */

                break;
              }

              case SUPP_LINE_BISECTOR: {
                /* extract two supporting lines */
                Rat_line_2 supp_line1; Rat_line_2 supp_line2;
                CGAL_assertion(CGAL::assign(supp_line1, o1));
                CGAL_assertion(CGAL::assign(supp_line2, o2));

                /* create bisector, orient it according to curr_direction, make
                 * it into a ray */

                /* find actual next intersection of ray */

                /* get segment */

                break;
              }

              case ENDPOINT_BISECTOR: {
                /* extract two endpoints */
                Rat_point_2 endpoint1; Rat_point_2 endpoint2;
                CGAL_assertion(CGAL::assign(endpoint1, o1));
                CGAL_assertion(CGAL::assign(endpoint2, o2));

                /* create bisector, orient it according to curr_direction, make
                 * it into a ray */

                /* find actual next intersection of ray */

                /* get segment */

                break;
              }

              default: break; // should never happen
            }
            //TODO finish

            /* add the piece of the bisector to the OutputIterator o, update the
             * curr_pt to be the next intersection found (corrected when
             * determining the actual correct piece of the bisector) */
            //TODO
            //TODO REMEMBER TO UPDATE curr_pt
            //TODO REMEMBER TO UPDATE curr_direction

            break; //TODO remove (to avoid infinite loop)
          }

          /* start from one point, find which segment is to consider as partial
           * directrix of parabolic arc. It is the segment that has as endpoints
           * the point whose orthogonal line intersects the ray to create the
           * ray start point. This information is saved in the info.
           * Also find the focus of the parabola supporting the parabolic arc,
           * which is the closest endpoint of the other segment. */
          // Rat_line_2 directrix;
          // Rat_point_2 focus;
          // Parabola supporting_conic;
          // if (start_ray_info.second == S1_SOURCE || start_ray_info.second == S1_TARGET) {
          //   directrix = s1.supporting_line();
          //   focus =
          //   (sqdistance(s2.source(), start_pt) < sqdistance(s2.target(), start_pt)) ?
          //   s2.source() : s2.target();
          //   supporting_conic = Parabola(directrix, focus);
          // }
          // else {
          //   directrix = s2.supporting_line();
          //   focus =
          //   (sqdistance(s1.source(), start_pt) < sqdistance(s1.target(), start_pt)) ?
          //   s1.source() : s1.target();
          //   supporting_conic = Parabola(directrix, focus);
          // }

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
         * In this case, the ray start points should be four, but only if the
         * intersection is not by one or two endpoints (weak intersection) */
        else {
          CGAL_assertion(CGAL::do_intersect(s1, s2)); // they HAVE to intersect
          CGAL_assertion(ray_info_list.size() == 4); //TODO correct for touching segments
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
