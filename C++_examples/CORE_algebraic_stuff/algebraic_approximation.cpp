#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <iostream>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

int main(int argc, char const *argv[]) {

  Algebraic a("1.9999999999999999999999999999999999999999999999999999999999999999999999");

  std::cout << "Created number a = " << a.toString() << '\n';
  std::cout << "Approximated with BigIntValue: " << a.BigIntValue() << '\n';
  return 0;
}

#endif
