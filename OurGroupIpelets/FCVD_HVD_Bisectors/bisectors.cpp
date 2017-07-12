#define CGAL_SDG_VERBOSE
#undef CGAL_SDG_VERBOSE

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Bisector_Linf.h>
#include <CGAL/Segment_Delaunay_graph_site_2.h>

#include <CGAL/Polychain_2.h>

#include <CGAL/Linf2D_voronoi_traits_2.h>
#include <CGAL/L2_voronoi_traits_2.h>
#include <CGAL/L2_coarse_HVD_traits_2.h>
#include <CGAL/L2_coarse_FCVD_traits_2.h>
#include <CGAL/envelope_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/clustergrabber.h>

namespace CGAL_bisectors{

  typedef CGAL::Cartesian<double>                           Kernel;
  typedef CGAL::Delaunay_triangulation_2<Kernel>            Delaunay;
  typedef CGAL::Regular_triangulation_euclidean_traits_2
    <Kernel,Kernel::FT>                                     RGt;
  typedef CGAL::Regular_triangulation_2<RGt>                Regular;
  //Apollonius
  typedef CGAL::Apollonius_graph_traits_2<Kernel>           AT;
  typedef CGAL::Apollonius_graph_2<AT>                      Apollonius;
  typedef Apollonius::Site_2                                ASite;
  // L_infinity bisector
  struct Gt_inf
    : public CGAL::Segment_Delaunay_graph_Linf_traits_2<
      Kernel,CGAL::Field_with_sqrt_tag> {};
  typedef Gt_inf::Site_2                                    Site_2;
  typedef
    CGAL::SegmentDelaunayGraphLinf_2::Bisector_Linf<Gt_inf> Inf_bis;
  // --------------------------------------------------------------------

  Inf_bis bisector_linf;

  // Voronoi diagrams
  typedef CGAL::Exact_predicates_exact_constructions_kernel VD_Kernel;
  typedef VD_Kernel::FT                                   Number_type;
  typedef VD_Kernel::Iso_rectangle_2                      Iso_rectangle_2;
  typedef VD_Kernel::Point_2                              VD_Point_2;
  typedef std::vector<VD_Point_2>                         Points;

  typedef CGAL::Linf2D_voronoi_traits_2<VD_Kernel>        VD_Traits_3;
  typedef VD_Traits_3::Surface_3   VD_Surface_3;
  typedef CGAL::Envelope_diagram_2<VD_Traits_3>           VD_Envelope_diagram_2;

  typedef CGAL::L2_voronoi_traits_2<VD_Kernel>   L2_VD_Traits_3;
  typedef L2_VD_Traits_3::Surface_3   L2_VD_Surface_3;
  typedef CGAL::Envelope_diagram_2<L2_VD_Traits_3> L2_VD_Envelope_diagram_2;

  typedef CGAL::L2_HVD_traits_2<VD_Kernel>       HVD_Traits_3;
  typedef HVD_Traits_3::Surface_3                HVD_Surface_3;
  typedef CGAL::Envelope_diagram_2<HVD_Traits_3> HVD_Envelope_diagram_2;

  typedef CGAL::L2_FCVD_traits_2<VD_Kernel>       FCVD_Traits_3;
  typedef FCVD_Traits_3::Surface_3                FCVD_Surface_3;
  typedef CGAL::Envelope_diagram_2<FCVD_Traits_3> FCVD_Envelope_diagram_2;

typedef  HVD_Envelope_diagram_2 Envelope;

  const unsigned int num_entries = 12;

  const std::string sublabel[] = {
    "two points euclidean bisector",
    "two points L_inf bisector",
    "point/segment L_inf-parabola",
    "two sites L_inf bisector",
    "Linf 2D Voronoi Diagram",
    "L2 farthest Voronoi Diagram",
    "Hausdorff Voronoi Diagram",
    "L2 FVD with polygonal input",
    "L2 farthest color Voronoi diagram (FCVD)",
    "L2 NVD with polygonal input",
    "L2 FCVDstar",
    "L2 farthest segment Voronoi Diagram",
    "Help"
  };

  const std::string helpmsg[] = {
    "Draw the L2 bisector of two points",
    "Draw the L_inf bisector of two points",
    "Draw the L_inf parabola for a point and a segment",
    "Draw the L_inf bisector for two sites (point or segment)",
    "Draw the L_inf Voronoi diagram for points with 2D bisectors",
    "Draw the L2 farthest Voronoi diagram for points",
    "Draw the Hausdorff Voronoi diagram for points",
    "Draw the L2 FVD for points of given clusters",
    "Draw the farthest color Voronoi diagram for points",
    "Draw the L2 NVD for points of given clusters",
    "Draw the L2 FCVDstar for points of given clusters",
    "Draw the L2 farthest Voronoi diagram for segments",
  };

  class bisectorIpelet
    : public CGAL::Ipelet_base<Kernel,num_entries> {
      public:
        bisectorIpelet()
          :CGAL::Ipelet_base<Kernel,num_entries>
             ("Bisectors",sublabel,helpmsg){}
        void protected_run(int);

        // cluster grabber

        template <class output_iterator>
        struct
        Cluster_grabber:public
        CGAL::internal::Cluster_grabber<Kernel,output_iterator>{
          Cluster_grabber(output_iterator it):
            CGAL::internal::Cluster_grabber<Kernel,output_iterator>(it){}
        };

        template<class output_iterator>
        boost::function_output_iterator<Cluster_grabber<output_iterator> >
        cluster_grabber(output_iterator it){
          return boost::make_function_output_iterator(
              Cluster_grabber<output_iterator>(it));
        }

    };
  // --------------------------------------------------------------------

  void bisectorIpelet::protected_run(int fn)
  {
    Delaunay dt;     //Voronoi of points
    Regular rt;     //power diagram
    Apollonius apo;     //apollonius

    if (fn == (num_entries-1)) {
      show_help();
      return;
    }

    std::list<Point_2> pt_list;
    std::list<Segment_2> sg_list;
    std::list<VD_Point_2> vd_pt_list;

    typedef CGAL::Polygon_2<Kernel>               Cluster_2;
    typedef CGAL::Polygon_2<VD_Kernel>            VD_Cluster_2;

    std::list<Cluster_2>    cluster_list;
    std::list<VD_Cluster_2> vd_cluster_list;

    Iso_rectangle_2 bbox;

    if ((fn == 6) or (fn == 7) or (fn == 8) or (fn == 9) or (fn == 10)) {
      // HVD, FVD from clusters, FCVD, NVD, FCVDstar from clusters
      // use cluster grabber:
      // a 1-point cluster {p}   is denoted by a point p (mark in ipe);
      // a 2-point cluster {p,q} is denoted by a segment pq;
      // a cluster with more than two points
      //   is denoted by a polygon with these points as vertices
      bbox =
        read_active_objects(
          CGAL::dispatch_or_drop_output
          <Point_2,Segment_2,Polygon_2>(
          cluster_grabber(std::back_inserter(cluster_list)),
          cluster_grabber(std::back_inserter(cluster_list)),
          cluster_grabber(std::back_inserter(cluster_list))
          )
          );
    } else { // all other cases
      bbox =
        read_active_objects(
          CGAL::dispatch_or_drop_output
          <Point_2,Polygon_2,Segment_2>(
          std::back_inserter(pt_list),
          segment_grabber(std::back_inserter(sg_list)),
          std::back_inserter(sg_list)
          )
          );
    }



    switch(fn){
      case 0:
        if (pt_list.empty()){
          print_error_message(("No mark selected"));
          return;
        }

        if (pt_list.size() != 2){
          print_error_message(("Exactly two points should be selected"));
          return;
        }

        break;
      case 1:
        // BISECTOR L_INF for two points
        if (pt_list.empty()){
          print_error_message(("No mark selected"));
          return;
        }

        if(pt_list.size() != 2) {
          print_error_message(("Exactly two points should be selected"));
          return;
        }

        break;

      case 2:
        // L_INF PARABOLA (one point and one segment)
        if (pt_list.empty() && sg_list.empty()){
          print_error_message(("No mark, no segment and no polygon selected"));
          return;
        }

        if(pt_list.size() + sg_list.size() != 2) {
          print_error_message(("Exactly two components should be selected"));
          return;
        }

        if(pt_list.size() != 1) {
          print_error_message(
              ("Exactly one point and one segment should be selected"));
          return;
        }

        break;

      case 3:
        // L_INF bisector of two sites (a site can be a point or a segment)
        if (pt_list.empty() && sg_list.empty()){
          print_error_message(("No mark, no segment and no polygon selected"));
          return;
        }

        if(pt_list.size() + sg_list.size() != 2) {
          print_error_message(("Exactly two components should be selected"));
          return;
        }

        break;

      case 4:
      case 5:
        {
          std::list<Point_2>::iterator it;
          for(it = pt_list.begin(); it != pt_list.end(); ++it) {
            vd_pt_list.push_back(VD_Point_2(it->x(), it->y()));
          }

          //Voronoi diagram for points.
          if (vd_pt_list.empty()){
            print_error_message(("No mark selected"));
            return;
          }
        }
        break;

      case 6:
      case 8:
      //  std::cout << "debug HVD/FCVD input check" << std::endl;
        if (cluster_list.empty()) {
          print_error_message(("No cluster selected"));
          return;
        } else {
          std::list<Cluster_2>::iterator clusterit;
          for (clusterit = cluster_list.begin();
               clusterit != cluster_list.end();
               clusterit++) {
            VD_Cluster_2 tempcluster;
            Cluster_2::Vertex_iterator it;
            for (it = clusterit->vertices_begin();
                 it != clusterit->vertices_end();
                 it++)
            {
              tempcluster.push_back(VD_Point_2(it->x(), it->y()));
            } // individual cluster for end
            vd_cluster_list.push_back(tempcluster);
          }
	  // family of clusters for end

        }
        break;

      case 7:
      case 9:
       // std::cout << "debug FVD/NVD cluster input check" << std::endl;
        if (cluster_list.empty()) {
          print_error_message(("No points selected"));
          return;
        } else {
          if (cluster_list.size() > 2) {
            print_error_message(("Too many polygons selected"));
            return;
          }

          std::list<Cluster_2>::iterator clusterit;
          for (clusterit = cluster_list.begin();
               clusterit != cluster_list.end();
               clusterit++) {
            Cluster_2::Vertex_iterator it;
            for (it = clusterit->vertices_begin();
                 it != clusterit->vertices_end();
                 it++)
            {
              vd_pt_list.push_back(VD_Point_2(it->x(), it->y()));
            }
          }

          //Voronoi diagram for points.
          if (vd_pt_list.empty()){
            print_error_message(("No mark selected"));
            return;
          }
        }

        break;
      case 10:
         if (cluster_list.empty()) {
                     print_error_message(("No points selected"));
		     }
         break;
	 }
     //end of switch

    Kernel::FT incr_len=(fn<2)?50:75;
    //slightly increase the size of the Bbox
    bbox = Iso_rectangle_2(
      bbox.min()+Kernel::Vector_2(-incr_len,-incr_len),
      bbox.max()+Kernel::Vector_2(incr_len,incr_len));

    if (fn == 0) {
      std::list<Point_2>::iterator it;
      it = pt_list.begin();
      Point_2 p = *it;
      ++it;
      Point_2 q = *it;
      Point_2 midp = CGAL::midpoint(p, q);
      Line_2 l(q,p);
      draw_in_ipe(l.perpendicular(midp));
    }

    if(fn == 1){
      std::list<Point_2>::iterator it;
      it = pt_list.begin();
      Point_2 p = *it;
      Site_2 sp = Site_2::construct_site_2(p);
      ++it;
      Point_2 q = *it;
      Site_2 sq = Site_2::construct_site_2(q);
      Inf_bis::Polychainline pcl = bisector_linf(sp, sq) ;
      Inf_bis::Polychainline::Vertex_const_iterator it1 = pcl.vertices_begin();
      if(pcl.size() == 1) {
        Point_2 firstpt = *it1;
        CGAL::Direction_2<Kernel> incomingDir = pcl.get_incoming();
        CGAL::Direction_2<Kernel> outgoingDir = pcl.get_outgoing();
        draw_in_ipe(Ray_2(firstpt, incomingDir));
        draw_in_ipe(Ray_2(firstpt, outgoingDir));

      }

      else if(pcl.size() == 2) {
        Point_2 firstpt = *it1;
        ++it1;
        Point_2 lastpt = *it1;
        draw_in_ipe(Ray_2(firstpt, pcl.get_incoming()));
        draw_in_ipe(Segment_2(firstpt, lastpt));
        draw_in_ipe(Ray_2(lastpt, pcl.get_outgoing()));

      }

      else {
        Point_2 firstpt = *it1;
        Ray_2 r1(firstpt, pcl.get_incoming());
        draw_in_ipe(r1);
        Inf_bis::Polychainline::Vertex_const_iterator it2 = it1+1;
        for(; it2!=pcl.vertices_end(); ++it1, ++it2) {
          draw_in_ipe(Segment_2(*it1, *it2));
        }
        draw_in_ipe(Ray_2(*it1, pcl.get_outgoing()));
      }

    }

    if(fn==2) {
      std::list<Point_2>::iterator it;
      it = pt_list.begin();
      Point_2 p = *it;
      Site_2 sp = Site_2::construct_site_2(p);

      std::list<Segment_2>::iterator it1;
      it1 = sg_list.begin();
      Segment_2 seg = *it1;
      Site_2 sq = Site_2::construct_site_2(seg.source(), seg.target());

      Inf_bis::Polychainline pcl = bisector_linf(sp, sq) ;
      Inf_bis::Polychainline::Vertex_const_iterator it2 = pcl.vertices_begin();
      Point_2 pt = *it2;
      draw_in_ipe(Ray_2(pt, pcl.get_incoming()));
      Inf_bis::Polychainline::Vertex_const_iterator it3 = it2+1;
      for(; it3!=pcl.vertices_end(); ++it2, ++it3) {
        draw_in_ipe(Segment_2(*it2, *it3));
      }
      draw_in_ipe(Ray_2(*it2, pcl.get_outgoing()));

    }

    if(fn==3) {
      std::list<Point_2>::iterator ptIt;
      ptIt = pt_list.begin();

      std::list<Segment_2>::iterator segIt;
      segIt = sg_list.begin();

      Site_2 sp, sq;
      //PP
      if(pt_list.size() <= 2 and sg_list.empty()) {
        Point_2 p = *ptIt;
        sp = Site_2::construct_site_2(p);
        ++ptIt;
        Point_2 q = *ptIt;
        sq = Site_2::construct_site_2(q);

      }
      //PS and SP
      else if (pt_list.size() == 1 and sg_list.size() == 1) {
        Point_2 p = *ptIt;
        sp = Site_2::construct_site_2(p);
        Segment_2 q = *segIt;
        sq = Site_2::construct_site_2(q.source(), q.target());
      }
    //SS
      else if (pt_list.empty() and sg_list.size() == 2) {
        Segment_2 p = *segIt;
        sp = Site_2::construct_site_2(p.source(), p.target());
        ++segIt;
        Segment_2 q = *segIt;
        sq = Site_2::construct_site_2(q.source(), q.target());
      }
        Inf_bis::Polychainline pcl = bisector_linf(sp, sq) ;
        Inf_bis::Polychainline::Vertex_const_iterator it2 = pcl.vertices_begin();
        Point_2 pt = *it2;
        draw_in_ipe(Ray_2(pt, pcl.get_incoming()));
        Inf_bis::Polychainline::Vertex_const_iterator it3 = it2+1;
        for(; it3!=pcl.vertices_end(); ++it2, ++it3) {
          draw_in_ipe(Segment_2(*it2, *it3));
        }
        draw_in_ipe(Ray_2(*it2, pcl.get_outgoing()));

    }

    if (fn==4) {
      VD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new VD_Envelope_diagram_2();
      CGAL::lower_envelope_3
        (vd_pt_list.begin(), vd_pt_list.end(), *m_envelope_diagram);

      //computes the bounding box
      VD_Point_2 bottom_left (bbox.min().x(), bbox.min().y());
      VD_Point_2 top_right (bbox.max().x(), bbox.max().y());

      for(VD_Envelope_diagram_2::Vertex_const_iterator vit =
            m_envelope_diagram->vertices_begin();
          vit != m_envelope_diagram->vertices_end();
          vit++) {
          VD_Point_2 vp = VD_Point_2(vit->point());
        if(CGAL::compare(vp.x(), bottom_left.x()) == CGAL::SMALLER) bottom_left = VD_Point_2(vp.x(), bottom_left.y());
        if(CGAL::compare(vp.y(), bottom_left.y()) == CGAL::SMALLER) bottom_left = VD_Point_2(bottom_left.x(), vp.y());
        if(CGAL::compare(vp.x(), top_right.x()) == CGAL::LARGER) top_right = VD_Point_2(vp.x(), top_right.y());
        if(CGAL::compare(vp.y(), top_right.y()) == CGAL::LARGER) top_right = VD_Point_2(top_right.x(), vp.y());

         //if one wants to display vertices of the VD as well, that's it
        //Point_2 p (to_double(vp.x()), to_double(vp.y()));
        //draw_in_ipe(p, bbox);
      }

      Point_2 bl (to_double(bottom_left.x()), to_double(bottom_left.y()));
      Point_2 tr (to_double(top_right.x()), to_double(top_right.y()));

      Kernel::FT incr_len= 50;

      bbox = Iso_rectangle_2(
                    bl + Kernel::Vector_2(-incr_len,-incr_len),
                    tr + Kernel::Vector_2(incr_len,incr_len));

      // draws edges
      for(VD_Envelope_diagram_2::Edge_const_iterator eit = m_envelope_diagram->edges_begin();
      eit != m_envelope_diagram->edges_end(); eit++) {
        if(eit->curve().is_segment()){
                Point_2 p1 (to_double(eit->curve().segment().source().x()), to_double(eit->curve().segment().source().y()));
                Point_2 p2 (to_double(eit->curve().segment().target().x()), to_double(eit->curve().segment().target().y()));
                draw_in_ipe(Segment_2(p1, p2), bbox);
        }else if (eit->curve().is_ray()){
                Point_2 p (to_double(eit->curve().ray().source().x()), to_double(eit->curve().ray().source().y()));
                CGAL::Direction_2<Kernel> d (to_double(eit->curve().ray().direction().dx()), to_double(eit->curve().ray().direction().dy()));
                draw_in_ipe(Ray_2(p, d), bbox);
        }else if(eit->curve().is_line()){
                Line_2 l (to_double(eit->curve().line().a()), to_double(eit->curve().line().b()), to_double(eit->curve().line().c()));
                draw_in_ipe(l, bbox);
        }
      }
    } // end of case: fn==4

    if ((fn==5) or (fn==7) or (fn==9)) {
      L2_VD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new L2_VD_Envelope_diagram_2();
      if ((fn==5) or (fn==7)) {
        CGAL::upper_envelope_3
          (vd_pt_list.begin(), vd_pt_list.end(), *m_envelope_diagram);
      } else {
        // NVD
        CGAL::lower_envelope_3
          (vd_pt_list.begin(), vd_pt_list.end(), *m_envelope_diagram);
      }

      //computes the bounding box
      VD_Point_2 bottom_left (bbox.min().x(), bbox.min().y());
      VD_Point_2 top_right (bbox.max().x(), bbox.max().y());

      for(L2_VD_Envelope_diagram_2::Vertex_const_iterator vit =
            m_envelope_diagram->vertices_begin();
          vit != m_envelope_diagram->vertices_end();
          vit++) {
          VD_Point_2 vp = VD_Point_2(vit->point());
        if(CGAL::compare(vp.x(), bottom_left.x()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(vp.x(), bottom_left.y());
        if(CGAL::compare(vp.y(), bottom_left.y()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(bottom_left.x(), vp.y());
        if(CGAL::compare(vp.x(), top_right.x()) == CGAL::LARGER)
          top_right = VD_Point_2(vp.x(), top_right.y());
        if(CGAL::compare(vp.y(), top_right.y()) == CGAL::LARGER)
          top_right = VD_Point_2(top_right.x(), vp.y());

        //if one wants to display vertices of the VD as well, that's it
        //Point_2 p (to_double(vp.x()), to_double(vp.y()));
        //draw_in_ipe(p, bbox);
      }

      Point_2 bl (to_double(bottom_left.x()), to_double(bottom_left.y()));
      Point_2 tr (to_double(top_right.x()), to_double(top_right.y()));

      Kernel::FT incr_len = 50;

      bbox = Iso_rectangle_2(
                    bl + Kernel::Vector_2(-incr_len,-incr_len),
                    tr + Kernel::Vector_2(incr_len,incr_len));

      // draws edges
      for(L2_VD_Envelope_diagram_2::Edge_const_iterator eit =
            m_envelope_diagram->edges_begin();
          eit != m_envelope_diagram->edges_end();
          eit++) {
        if (eit->curve().is_segment()) {
          Point_2 p1 (to_double(eit->curve().segment().source().x()),
                      to_double(eit->curve().segment().source().y()));
          Point_2 p2 (to_double(eit->curve().segment().target().x()),
                      to_double(eit->curve().segment().target().y()));
          draw_in_ipe(Segment_2(p1, p2), bbox);
        } else if (eit->curve().is_ray()) {
          Point_2 p (to_double(eit->curve().ray().source().x()),
                     to_double(eit->curve().ray().source().y()));
          CGAL::Direction_2<Kernel> d
            (to_double(eit->curve().ray().direction().dx()),
             to_double(eit->curve().ray().direction().dy()));
          draw_in_ipe(Ray_2(p, d), bbox);
        } else if (eit->curve().is_line()) {
          Line_2 l (to_double(eit->curve().line().a()),
                    to_double(eit->curve().line().b()),
                    to_double(eit->curve().line().c()));
          draw_in_ipe(l, bbox);
        }
      } // end of draw edges
    } // end of case: fn==5

    if (fn==6) { // HVD
      HVD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new HVD_Envelope_diagram_2();


      CGAL::lower_envelope_3
        (vd_cluster_list.begin(), vd_cluster_list.end(), *m_envelope_diagram);

      //computes the bounding box
      VD_Point_2 bottom_left (bbox.min().x(), bbox.min().y());
      VD_Point_2 top_right (bbox.max().x(), bbox.max().y());

      for(HVD_Envelope_diagram_2::Vertex_const_iterator vit =
            m_envelope_diagram->vertices_begin();
          vit != m_envelope_diagram->vertices_end();
          vit++) {
          VD_Point_2 vp = VD_Point_2(vit->point());
        if(CGAL::compare(vp.x(), bottom_left.x()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(vp.x(), bottom_left.y());
        if(CGAL::compare(vp.y(), bottom_left.y()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(bottom_left.x(), vp.y());
        if(CGAL::compare(vp.x(), top_right.x()) == CGAL::LARGER)
          top_right = VD_Point_2(vp.x(), top_right.y());
        if(CGAL::compare(vp.y(), top_right.y()) == CGAL::LARGER)
          top_right = VD_Point_2(top_right.x(), vp.y());

        //if one wants to display vertices of the VD as well, that's it
        //Point_2 p (to_double(vp.x()), to_double(vp.y()));
        //draw_in_ipe(p, bbox);
      }

      Point_2 bl (to_double(bottom_left.x()), to_double(bottom_left.y()));
      Point_2 tr (to_double(top_right.x()), to_double(top_right.y()));

      Kernel::FT incr_len = 50;

      bbox = Iso_rectangle_2(
                    bl + Kernel::Vector_2(-incr_len,-incr_len),
                    tr + Kernel::Vector_2(incr_len,incr_len));

      // draws edges
      for(HVD_Envelope_diagram_2::Edge_const_iterator eit =
            m_envelope_diagram->edges_begin();
          eit != m_envelope_diagram->edges_end();
          eit++) {
        if (eit->curve().is_segment()) {
          Point_2 p1 (to_double(eit->curve().segment().source().x()),
                      to_double(eit->curve().segment().source().y()));
          Point_2 p2 (to_double(eit->curve().segment().target().x()),
                      to_double(eit->curve().segment().target().y()));
          draw_in_ipe(Segment_2(p1, p2), bbox);
        } else if (eit->curve().is_ray()) {
          Point_2 p (to_double(eit->curve().ray().source().x()),
                     to_double(eit->curve().ray().source().y()));
          CGAL::Direction_2<Kernel> d
            (to_double(eit->curve().ray().direction().dx()),
             to_double(eit->curve().ray().direction().dy()));
          draw_in_ipe(Ray_2(p, d), bbox);
        } else if (eit->curve().is_line()) {
          Line_2 l (to_double(eit->curve().line().a()),
                    to_double(eit->curve().line().b()),
                    to_double(eit->curve().line().c()));
          draw_in_ipe(l, bbox);
        }
      } // end of draw edges
    } // end of case: fn==6

    if (fn==8) { // FCVD
      FCVD_Envelope_diagram_2 *m_envelope_diagram;
      m_envelope_diagram = new FCVD_Envelope_diagram_2();

      CGAL::upper_envelope_3
        (vd_cluster_list.begin(), vd_cluster_list.end(), *m_envelope_diagram);

      //computes the bounding box
      VD_Point_2 bottom_left (bbox.min().x(), bbox.min().y());
      VD_Point_2 top_right (bbox.max().x(), bbox.max().y());

      for(FCVD_Envelope_diagram_2::Vertex_const_iterator vit =
            m_envelope_diagram->vertices_begin();
          vit != m_envelope_diagram->vertices_end();
          vit++) {
          VD_Point_2 vp = VD_Point_2(vit->point());
        if(CGAL::compare(vp.x(), bottom_left.x()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(vp.x(), bottom_left.y());
        if(CGAL::compare(vp.y(), bottom_left.y()) == CGAL::SMALLER)
          bottom_left = VD_Point_2(bottom_left.x(), vp.y());
        if(CGAL::compare(vp.x(), top_right.x()) == CGAL::LARGER)
          top_right = VD_Point_2(vp.x(), top_right.y());
        if(CGAL::compare(vp.y(), top_right.y()) == CGAL::LARGER)
          top_right = VD_Point_2(top_right.x(), vp.y());

        //if one wants to display vertices of the VD as well, that's it
        //Point_2 p (to_double(vp.x()), to_double(vp.y()));
        //draw_in_ipe(p, bbox);
      }

      Point_2 bl (to_double(bottom_left.x()), to_double(bottom_left.y()));
      Point_2 tr (to_double(top_right.x()), to_double(top_right.y()));

      Kernel::FT incr_len = 50;

      bbox = Iso_rectangle_2(
                    bl + Kernel::Vector_2(-incr_len,-incr_len),
                    tr + Kernel::Vector_2(incr_len,incr_len));

      // draws edges
      for(FCVD_Envelope_diagram_2::Edge_const_iterator eit =
            m_envelope_diagram->edges_begin();
          eit != m_envelope_diagram->edges_end();
          eit++) {
        if (eit->curve().is_segment()) {
          Point_2 p1 (to_double(eit->curve().segment().source().x()),
                      to_double(eit->curve().segment().source().y()));
          Point_2 p2 (to_double(eit->curve().segment().target().x()),
                      to_double(eit->curve().segment().target().y()));
          draw_in_ipe(Segment_2(p1, p2), bbox);
        } else if (eit->curve().is_ray()) {
          Point_2 p (to_double(eit->curve().ray().source().x()),
                     to_double(eit->curve().ray().source().y()));
          CGAL::Direction_2<Kernel> d
            (to_double(eit->curve().ray().direction().dx()),
             to_double(eit->curve().ray().direction().dy()));
          draw_in_ipe(Ray_2(p, d), bbox);
        } else if (eit->curve().is_line()) {
          Line_2 l (to_double(eit->curve().line().a()),
                    to_double(eit->curve().line().b()),
                    to_double(eit->curve().line().c()));
          draw_in_ipe(l, bbox);
        }
      } // end of draw edges
    } // end of case: fn==8
  }
  // end of void bisectorIpelet::protected_run(int fn)

}

CGAL_IPELET(CGAL_bisectors::bisectorIpelet)
