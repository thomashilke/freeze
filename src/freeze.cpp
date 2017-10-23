#include <iostream>

#include "freeze.hpp"

int main(int argc, char *argv[]) {
  try {
    if (argc != 2)
      std::cout << "Usage: " << argv[0] << "  <configuration file>" << std::endl;
    else {
      parameter::collection p;
      p.read_from_file(argv[1]);

      for (std::size_t i(0); i < p.get_collection_size(); ++i) {
        std::cout << "running parameter set collection #" << i + 1 << " of " << p.get_collection_size() << std::endl;        
        p.set_current_collection(i);
        
        freeze f(p);
        f.run();
      }
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
