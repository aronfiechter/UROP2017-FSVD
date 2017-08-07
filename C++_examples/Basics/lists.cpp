#include <iostream>
#include <list>
#include <iterator>

/* function that inserts */
void list_inserter(std::back_insert_iterator<std::list<int>> result);
void list_inserter(std::back_insert_iterator<std::list<int>> result) {
  *result++ = 1;
  *result++ = 2;
  *result++ = 3;
  *result++ = 4;
  *result++ = 5;
}

int main() {

  std::list<int> l;

  list_inserter(std::back_insert_iterator<std::list<int>>(l));

  std::cout << "Values in l:";
  for (std::list<int>::iterator it = l.begin(); it != l.end(); ++it) {
    std::cout << ' ' << *it;
  }
  std::cout << std::endl;

}
