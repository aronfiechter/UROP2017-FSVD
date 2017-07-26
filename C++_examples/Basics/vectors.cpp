#include <iostream>
#include <vector>

int main() {

  std::vector<int> v;

  v.push_back(1);
  v.push_back(2);
  v.push_back(3);
  v.push_back(4);
  v.push_back(5);

  std::cout << "Values in v:";
  for (std::vector<int>::iterator it = v.begin(); it != v.end(); ++it) {
    std::cout << ' ' << *it;
  }
  std::cout << '\n';

}
