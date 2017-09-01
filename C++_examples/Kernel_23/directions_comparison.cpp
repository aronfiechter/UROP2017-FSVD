#include <iostream>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Line_2 Line_2;
typedef Kernel::Direction_2 Direction_2;

int main() {
  Point_2 p(0, 0), q(1, 1), r(1, 0);

  std::cout << "p = " << p << std::endl;
  std::cout << "q = " << q.x() << " " << q.y() << std::endl;

  Line_2 l(p, q);
  std::cout << "l = " << l << " with direction " << l.direction() << std::endl;

  Line_2 flat(p, r);
  std::cout << "flat = " << flat << " with direction " << flat.direction() << std::endl;

  Line_2 bisector = CGAL::bisector(l, flat);
  std::cout << "bisector(l, flat) = " << bisector << " with direction " << bisector.direction() << std::endl;

  Line_2 bisector_2 = CGAL::bisector(l.opposite(), flat);
  std::cout << "bisector(-l, flat) = " << bisector_2 << " with direction " << bisector_2.direction() << std::endl;

  return 0;
}
