#include <iostream>
#include <vector>

void print_i(int i);

int main() {

  std::vector<int> v{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  for (auto& x : v) {
    ++x;
  }

  std::for_each(v.begin(), v.end(), print_i);

  std::cout << std::endl;

  return 0;
}

void print_i(int i) {
  std::cout << i << ' ';
}
