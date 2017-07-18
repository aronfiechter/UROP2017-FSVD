// created by Aron Fiechter on 2017-07-13.
// This file is a copy of L2_voronoi_traits_2, except that it is for segments
// instead of points.
// We import CGAL/Parabola_segment_2.h because we need it for segment bisectors.

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

#include <CGAL/Parabola_segment_2.h>

namespace CGAL {

template <class Kernel_>
class L2_segment_voronoi_traits_2 : public Arr_linear_traits_2<Kernel_> {

public:
  typedef Kernel_                                           Kernel;
  typedef Arr_linear_traits_2<Kernel>                       Base;
  typedef typename Kernel::FT                               FT;
  typedef typename Base::Multiplicity                       Multiplicity;

  typedef typename Base::Point_2                            Point_2;
  typedef typename Base::Curve_2                            Curve_2;
  typedef typename Base::X_monotone_curve_2                 X_monotone_curve_2;
  typedef typename Kernel::Segment_2                        Segment_2;
  typedef typename Kernel::Ray_2                            Ray_2;
  typedef typename Kernel::Line_2                           Line_2;
  typedef typename Kernel::Direction_2                      Direction_2;

  /* For parts of segment bisectors */
  typedef CGAL::Parabola_segment_2<Kernel>                  Parabola_segment_2;

  /* Define segments as surfaces for simplicity. Each segment should be thought
   * of as a surface representing the function distance, with higher values on
   * the z axis mean farthest points */
  typedef Segment_2       Xy_monotone_surface_3;
  typedef Segment_2       Surface_3;



protected:
  typedef std::pair<X_monotone_curve_2, Multiplicity>       Intersection_curve;

  /* Returns the squared distance between two points in L2 metric. */
  static FT sqdistance(const Point_2& p1, const Point_2& p2) {
    FT diffx = p1.x() - p2.x();
    FT diffy = p1.y() - p2.y();
    FT sqdist = CGAL::square(diffx) + CGAL::square(diffy);
    return sqdist;
  }

  /* Returns the squared distance between a point and a segment in L2 metric. */
  static FT sqdistance(const Point_2& p, const Segment_2& s){
    /* find projection of p on supporting line of s */
    Line_2 l = s.supporting_line();
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
      /* otherwise, for now just make bisector of the two source points of the
       * segments, to test if this works */
      else {
        *o++ = CGAL::make_object(
          Intersection_curve(CGAL::bisector(s1.source(), s2.source()), 0)
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
      return CGAL::compare(sqdistance(p, h1), sqdistance(p, h2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::EQUAL;
    }

    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      /* if the two unbounded surfaces do not intersect, then they must
       * represent the same segment's distance function */
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

}; // class L2_segment_voronoi_traits_2
} // namespace CGAL

#endif // CGAL_L2_SEGMENT_VORONOI_TRAITS_2_H
