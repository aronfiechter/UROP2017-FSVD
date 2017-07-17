#include <iostream>

class Person {
public:
  Person(std::string _name, int _age, float _height)
    : name(_name), age(_age), height(_height) {}

  std::string get_name() {
    return name;
  }
  float get_height() {
    return height;
  }
  int get_age() {
    return age;
  }

private:
  std::string name;
  int age;
  float height;
};

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

  std::cout << "You created a person called " << me.get_name()
            << " of age " << me.get_age()
            << " and height " << me.get_height()
            << "." << std::endl;
}
