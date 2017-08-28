#include <iostream>
#include <list>

class Person {
public:
  Person(std::string _name, int _age, float _height)
    : name(_name), age(_age), height(_height) {}

  /* getters */
  std::string get_name(void) const {
    return name;
  }
  float get_height(void) const {
    return height;
  }
  int get_age(void) const {
    return age;
  }

  /* setters */
  void set_height(const float _height) {
    this->height = _height;
  }

private:
  std::string name;
  int age;
  float height;
};

std::ostream& operator<<(std::ostream& os, const Person& person) {
  os << "Name: " << person.get_name() << "\n"
     << "Age: " << person.get_age() << "\n"
     << "Height: " << person.get_height() << "\n"
  ;
  return os;
}

int main() {
  std::string my_name;
  int my_age;
  float my_height;
  std::cout << "Enter a name: ";
  std::cin >> my_name;
  std::cout << "Enter an age: ";
  std::cin >> my_age;
  std::cout << "Enter a height: ";
  std::cin >> my_height;
  Person me(my_name, my_age, my_height);
  Person mom("Mom", 56, 156);
  Person dad("Dad", 54, 178);

  std::cout << me << "\n" << mom << "\n" << dad << "\n";

  std::list<Person> people;
  people.push_back(me);
  people.push_back(mom);
  people.push_back(dad);

  dad.set_height(177);
  people.pop_back();
  people.push_back(dad);

  std::cout << "Updated height of dad:\n" << dad << '\n';

}
