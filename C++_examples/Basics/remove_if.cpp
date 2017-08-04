// remove_if example
#include <iostream>     // std::cout
#include <vector>
#include <algorithm>    // std::remove_if

bool IsOdd (int i);
bool IsOdd (int i) {
  return ((i % 2) == 1);
}

int main () {
  std::vector<int> myints = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  typename std::vector<int>::iterator new_end
    = std::remove_if (myints.begin(), myints.end(), IsOdd);

  myints.resize(static_cast<unsigned long>(new_end - myints.begin()));

  std::cout << "the range contains:";
  for (auto& x : myints) {
    std::cout << " " << x;
  }
  std::cout << '\n';

  return 0;
}
