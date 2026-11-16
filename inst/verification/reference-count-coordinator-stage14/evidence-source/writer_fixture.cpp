#include "preservation_core.h"
#include <iostream>
int main(int argc, char** argv) {
  if (argc != 6) return 2;
  try {
    mss::prepare(argv[1], argv[2], argv[3], std::stoull(argv[4]), std::stoi(argv[5]));
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
