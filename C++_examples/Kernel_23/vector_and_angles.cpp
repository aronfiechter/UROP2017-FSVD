#include <iostream>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_2 Line_2;
typedef Kernel::Vector_2 Vector_2;

int main()
{
  Vector_2 a(-3.1, 2.6), d(3.9, 0.9);

  std::cout << "a = " << a << ", ||a||^2 = " << a.squared_length() << std::endl;
  std::cout << "d = " << d << ", ||d||^2 = " << d.squared_length() << std::endl;

  std::cout << "Normalize..." << '\n';
  a = a / CGAL::sqrt(a.squared_length());
  d = d / CGAL::sqrt(d.squared_length());
  std::cout << "a = " << a << ", ||a||^2 = " << a.squared_length() << std::endl;
  std::cout << "d = " << d << ", ||d||^2 = " << d.squared_length() << std::endl;

  std::cout << "Get halfway vector..." << '\n';
  Vector_2 e = a + d;
  std::cout << "e = " << e << ", ||e||^2 = " << e.squared_length() << std::endl;

  std::cout << "Normalize..." << '\n';
  e = e / CGAL::sqrt(e.squared_length());
  std::cout << "e = " << e << ", ||e||^2 = " << e.squared_length() << std::endl;

  std::cout << "Get length of 4..." << '\n';
  e = 4 * e;
  std::cout << "e = " << e << ", ||e||^2 = " << e.squared_length() << std::endl;

  return 0;
}
