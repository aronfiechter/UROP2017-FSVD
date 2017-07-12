#include <iostream>
#include <cstring>

int main() {

  char s[1000];

  std::cin >> s;

  int Uppercase = 0; //modified by the lambda

  std::for_each(s, s + strlen(s), [&Uppercase] (char c) {
    if (isupper(c))
      Uppercase++;
  });

  std::cout << Uppercase << " uppercase letters in: " << s << std::endl;
}
