#include <iostream>
#include <string>

#include <cstring>

#include "cardnameC.h"

extern "C" {
void necard_();
}

int main(int argc, char const *argv[]) {

  if (argc < 2) {
    std::cout << "[ERROR]: Expected to be passed a card file." << std::endl;
    return 1;
  }

  std::strcpy(necardname_.fcard, argv[1]);

  std::string out_file = "neutescatvect.root";
  if (argc > 2) {
    out_file = argv[2];
  }

  necard_();
}
