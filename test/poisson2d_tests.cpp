#include "catch.hpp"
#include "../src/poisson2d.hpp"
#include "../src/jsl.hpp"

TEST_CASE("Case 5 - ADAD, GEnxn", "[Poisson2D]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(5);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  REQUIRE(p2d.N == 5);
  REQUIRE(p2d.x.size() == 25);
  REQUIRE(p2d.y.size() == 25);
  REQUIRE(p2d.u.size() == 25);

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  REQUIRE(M.size() == 15);
  REQUIRE(map.size() == 15);
  REQUIRE(b.size() == 15);

  std::vector<double> x(15), z(15);
  std::vector<size_t> ip(15);

  jsl::GEnxn<size_t, double, std::vector>(15, M, x, b, z, ip);
  jsl::map_vector<size_t, double, std::vector>(x, p2d.u, map);
  for (size_t i = 0; i < 5; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
    CHECK(p2d(i, 3) == Approx(p2d.x[i]));
    CHECK(p2d(i, 4) == Approx(p2d.x[i]));
  }

  // Test other stuff
}
