#include <iostream>

#include "freeze.hpp"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2)
      std::cout << "Usage: " << argv[0] << "  <configuration file>" << std::endl;
    else {
      parameter::collection p;
      p.read_from_file(argv[1]);
      
      freeze f(p);
      f.run();
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
