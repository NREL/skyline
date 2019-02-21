#include <vector>
#include <stdio.h>
#include "../include/skyline.hpp"

int main()
{
  // Solve the one-dimensional Laplace problem on the unit interval
  
  // Create the object and fill in the system of equations
  std::vector<size_t> heights{ {0, 1, 1, 1, 1, 1, 1} };
  skyline::SymmetricMatrix<size_t, float, std::vector> sky1(heights);
  for (size_t i = 0; i < 7; ++i) {
    sky1.diagonal(i) = 2.0;
  }
  for (size_t i = 0; i < 6; ++i) {
    sky1(i) = -1.0;
  }

  // Build the RHS
  std::vector<float> x{ {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0} };

  // Solve
  sky1.ldlt_solve(x);

  // Report out
  puts(" i      x          f");
  puts("--- ---------  ---------");
  for (size_t i = 0; i < 7; ++i) {
    printf("%2d  %.3e  %.3e\n", (int)(i+1), 0.125*(i+1), x[i]);
  }
  puts("--- ---------  ---------");

  exit(EXIT_SUCCESS);

}
