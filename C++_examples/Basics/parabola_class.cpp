#include <iostream>
#include <cassert>

class Parabola {
private: // coefficients of equation: rx^2 + sy^2 + txy + ux + vy + w = 0
  typedef int RT;
  RT _r; RT _s; RT _t; RT _u; RT _v; RT _w;

public:
  /* Construct using equation coefficients */
  Parabola(RT __r, RT __s, RT __t, RT __u, RT __v, RT __w)
    : _r(__r), _s(__s), _t(__t), _u(__u), _v(__v), _w(__w) {

    assert(_t * _t - 4 * _r * _s == 0); // curve is a parabola
  }

  /* Getters */
  RT r() { return _r; }
  RT s() { return _s; }
  RT t() { return _t; }
  RT u() { return _u; }
  RT v() { return _v; }
  RT w() { return _w; }
};

int main() {
  int r, s, t, u, v, w;
  std::cout << "Enter coefficients r s t u v w: ";
  std::cin >> r >> s >> t >> u >> v >> w;
  Parabola p(r, s, t, u, v, w);

  std::cout << std::endl << "Constructed parabola p: ";
  if (p.r() != 0) std::cout << p.r() << "x^2 + ";
  if (p.s() != 0) std::cout << p.s() << "y^2 + ";
  if (p.t() != 0) std::cout << p.t() << "xy + ";
  if (p.u() != 0) std::cout << p.u() << "x + ";
  if (p.v() != 0) std::cout << p.v() << "y + ";
  if (p.w() != 0) std::cout << p.w();
  std::cout << " = 0" << std::endl;
}
