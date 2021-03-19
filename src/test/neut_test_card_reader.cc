#include <fstream>
#include <iostream>

#include "CommonBlockIFace.h"

int main(int argc, char const *argv[]) {
  if (argc < 1) {
    std::cout << "[ERROR]: Expected the path of a card file to read."
              << std::endl;
  }
  neut::CommonBlockIFace::Initialize(argv[1]);

  if (argc > 2) {
    std::cout << "[INFO]: Writing common blocks to " << argv[2] << std::endl;
    std::ofstream fout(argv[2]);
    fout << neut::CommonBlockIFace::Get().ParamsToString() << std::endl;
  } else {
    std::cout << neut::CommonBlockIFace::Get().ParamsToString() << std::endl;
  }
}