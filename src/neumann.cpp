#include <iostream>

#include "freeze.hpp"


double sqr(double x) { return x * x; }

int main(int argc, char *argv[]) {
  try {
    if (argc != 2)
      std::cout << "Usage: " << argv[0] << "  <configuration file>" << std::endl;
    else {
      parameter::collection p;
      p.read_from_file(argv[1]);
      
      freeze f(p);
      auto solution(f.run());


      const double time_end(p.get_value<double>("time-end"));
      const double L2_error(
        std::sqrt(
          integrate<quad::edge::gauss5>(
            compose(sqr,
              make_expr(std::bind(neumann_exact_solution,
                                  time_end,
                                  std::placeholders::_1))
              - make_expr<freeze::fe_type>(solution)),
            solution.get_mesh())));

      std::cout << "L2_error(t = " << time_end << ") = " << L2_error << std::endl;
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
