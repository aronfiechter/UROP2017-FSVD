#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char const *argv[]) {

  char s[1000];

  cin >> s;

  int Uppercase = 0; //modified by the lambda

  for_each(s, s + strlen(s), [&Uppercase] (char c) {
    if (isupper(c))
      Uppercase++;
  });

  cout << Uppercase << " uppercase letters in: " << s << endl;
}
