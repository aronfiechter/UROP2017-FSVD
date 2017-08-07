#include <iostream>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Line_2 Line_2;

int main() {
  Point_2 p(0, 0), q(1, 1);

  std::cout << "p = " << p << std::endl;
  std::cout << "q = " << q.x() << " " << q.y() << std::endl;

  Line_2 l(p, q);
  std::cout << "l = " << l << " with direction " << l.direction() << std::endl;

  Line_2 rev_l = l.opposite();
  std::cout << "rev_l = " << rev_l << " with direction " << rev_l.direction() << std::endl;

  Point_2 m(1, 10);
  std::cout << "m = " << m << std::endl;

  std::cout << "m is on ";

  switch (l.oriented_side(m)){
  case CGAL::ON_ORIENTED_BOUNDARY:
    std::cout << "the line l" << std::endl;
    break;
  case CGAL::ON_POSITIVE_SIDE:
    std::cout << "the positive side of line l" << std::endl;
    break;
  case CGAL::ON_NEGATIVE_SIDE:
    std::cout << "the negative side of line l" << std::endl;
    break;
  }

  std::cout << "m is on ";

  switch (rev_l.oriented_side(m)){
  case CGAL::ON_ORIENTED_BOUNDARY:
    std::cout << "the line rev_l" << std::endl;
    break;
  case CGAL::ON_POSITIVE_SIDE:
    std::cout << "the positive side of line rev_l" << std::endl;
    break;
  case CGAL::ON_NEGATIVE_SIDE:
    std::cout << "the negative side of line rev_l" << std::endl;
    break;
  }

  return 0;
}
