// #define CGAL_SDG_VERBOSE
// #undef CGAL_SDG_VERBOSE
//
// #ifndef CGAL_SDG_VERBOSE
// #define CGAL_SDG_DEBUG(a)
// #else
// #define CGAL_SDG_DEBUG(a) { a }
// #endif

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

namespace CGAL_min_circle {

  // Voronoi diagrams
  typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  typedef CGAL::Min_circle_2_traits_2<K>                    Traits;
  typedef CGAL::Min_circle_2<Traits>                        Min_circle;
  typedef K::Point_2                                        Point;

  const unsigned int num_entries = 1;

  const std::string sublabel[] = {
    "Min circle enclosing points",
    "Help"
  };

  const std::string helpmsg[] = {
    "Draw the smallest circle enclosing the given points",
  };

  class minCircleIpelet
  : public CGAL::Ipelet_base<Kernel, num_entries> {
    public:
      // declare an ipelet called Enclosing circle, with 2 functions (including help message).
      minCircleIpelet()
        :CGAL::Ipelet_base<Kernel,num_entries>("Enclosing circle", sublabel, helpmsg) {}
      void protected_run(int);
  };


  void minCircleIpelet::protected_run(int fn) {

    switch (fn) {
    case 1:
      show_help(); // print an help message
      return;

    case 0:
      std::list<Point_2> pt_lst;

      // Recovering points using output iterator of type
      // Dispatch_or_drop_output_iterator
      read_active_objects(
        CGAL::dispatch_or_drop_output<Point_2>(std::back_inserter(pt_lst))
      );

      if (pt_lst.empty()) {
        print_error_message("No mark selected");
        return;
      }

      Delaunay dt;
      dt.insert(pt_lst.begin(),pt_lst.end());

      //draw the triangulation.
      draw_in_ipe(dt);
  };

  }

}

CGAL_IPELET(CGAL_bisectors::bisectorIpelet)
