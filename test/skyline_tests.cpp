#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../src/poisson2d.hpp"
#include "../src/skyline.hpp"
#include "../src/jsl.hpp"
#include <iostream>

TEST_CASE("Case 5 - ADAD, Skyline Incremental 4x3", "[Skyline]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(4,3);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  //for (size_t i = 0; i < M.size(); ++i) {
  //  for (size_t j = 0; j < M.size(); ++j) {
  //    std::cout << std::setw(2) << M[i][j] << ' ';
  //  }
  //  std::cout << "| " << b[i] << std::endl;
  //}

  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(M);

  REQUIRE(skyline.heights().size() == 6);
  CHECK(skyline.heights()[0] == 0);
  CHECK(skyline.heights()[1] == 0);
  CHECK(skyline.heights()[2] == 2);
  CHECK(skyline.heights()[3] == 2);
  CHECK(skyline.heights()[4] == 2);
  CHECK(skyline.heights()[5] == 2);

  REQUIRE(skyline.offsets().size() == 6);
  CHECK(skyline.offsets()[0] == 0);
  CHECK(skyline.offsets()[1] == 0);
  CHECK(skyline.offsets()[2] == 0);
  CHECK(skyline.offsets()[3] == 2);
  CHECK(skyline.offsets()[4] == 4);
  CHECK(skyline.offsets()[5] == 6);


  std::vector<double> au{ {-1.0, 0.0, -1.0, -1.0, -1.0, 0.0, -1.0, 0.0} };

  REQUIRE(skyline.upper().size() == 8);
  for (size_t i = 0; i < 8; ++i) {
    CHECK(skyline.upper()[i] == au[i]);
  }

  std::vector<double> v(6);
  jsl::UTDU<size_t, double, std::vector>(6, M, v);
  /*
  for (size_t i = 0; i < M.size(); ++i) {
    for (size_t j = 0; j < M.size(); ++j) {
      std::cout << std::setw(5) << M[i][j] << ' ';
    }
    std::cout << "| " << b[i] << std::endl;
  }*/

  skyline.utdu();
  for (size_t i = 0; i < 6; i++) {
    INFO("The index is " << i);
    CHECK(skyline.diagonal()[i] == M[i][i]);
  }
  skyline.forward_substitution(b);
  skyline.back_substitution(b);
  jsl::map_vector<size_t, double, std::vector>(b, p2d.u, map);
  for (size_t i = 0; i < 4; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
  }

}

TEST_CASE("Case 5 - ADAD, Skyline Incremental 4x4", "[Skyline]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(4, 4);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(M);

  REQUIRE(skyline.heights().size() == 8);
  CHECK(skyline.heights()[0] == 0);
  CHECK(skyline.heights()[1] == 0);
  CHECK(skyline.heights()[2] == 2);
  CHECK(skyline.heights()[3] == 2);
  CHECK(skyline.heights()[4] == 2);
  CHECK(skyline.heights()[5] == 2);
  CHECK(skyline.heights()[6] == 2);
  CHECK(skyline.heights()[7] == 2);

  REQUIRE(skyline.offsets().size() == 8);
  CHECK(skyline.offsets()[0] == 0);
  CHECK(skyline.offsets()[1] == 0);
  CHECK(skyline.offsets()[2] == 0);
  CHECK(skyline.offsets()[3] == 2);
  CHECK(skyline.offsets()[4] == 4);
  CHECK(skyline.offsets()[5] == 6);
  CHECK(skyline.offsets()[6] == 8);
  CHECK(skyline.offsets()[7] == 10);

  std::vector<double> au{ {-1.0, 0.0, -1.0, -1.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0, -1.0, 0.0} };

  REQUIRE(skyline.upper().size() == 12);
  for (size_t i = 0; i < 12; ++i) {
    INFO("The index is " << i);
    CHECK(skyline.upper()[i] == au[i]);
  }

  std::vector<double> v(8);
  jsl::UTDU<size_t, double, std::vector>(8, M, v);

  skyline.utdu();
  for (size_t i = 0; i < 6; i++) {
    INFO("The index is " << i);
    CHECK(skyline.diagonal()[i] == M[i][i]);
  }
  skyline.forward_substitution(b);
  skyline.back_substitution(b);
  jsl::map_vector<size_t, double, std::vector>(b, p2d.u, map);
  for (size_t i = 0; i < 4; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
    CHECK(p2d(i, 3) == Approx(p2d.x[i]));
  }

}

TEST_CASE("Case 5 - ADAD, Skyline Incremental 5x4", "[Skyline]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(5, 4);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(M);

  REQUIRE(skyline.heights().size() == 12);
  CHECK(skyline.heights()[0] == 0);
  CHECK(skyline.heights()[1] == 0);
  CHECK(skyline.heights()[2] == 0);
  CHECK(skyline.heights()[3] == 3);
  CHECK(skyline.heights()[4] == 3);
  CHECK(skyline.heights()[5] == 3);
  CHECK(skyline.heights()[6] == 3);
  CHECK(skyline.heights()[7] == 3);
  CHECK(skyline.heights()[8] == 3);
  CHECK(skyline.heights()[9] == 3);
  CHECK(skyline.heights()[10] == 3);
  CHECK(skyline.heights()[11] == 3);

  REQUIRE(skyline.offsets().size() == 12);
  CHECK(skyline.offsets()[0] == 0);
  CHECK(skyline.offsets()[1] == 0);
  CHECK(skyline.offsets()[2] == 0);
  CHECK(skyline.offsets()[3] == 0);
  CHECK(skyline.offsets()[4] == 3);
  CHECK(skyline.offsets()[5] == 6);

  std::vector<double> au{ {-1.0, 0.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0,
    -1.0, 0.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0 , -1.0, 0.0, 0.0 } };

  REQUIRE(skyline.upper().size() == 27);
  for (size_t i = 0; i < 8; ++i) {
    CHECK(skyline.upper()[i] == au[i]);
  }

  std::vector<double> v(12);
  jsl::UTDU<size_t, double, std::vector>(12, M, v);
  
  skyline.utdu();
  for (size_t i = 0; i < 12; i++) {
    INFO("The index is " << i);
    CHECK(skyline.diagonal()[i] == M[i][i]);
  }
  skyline.forward_substitution(b);
  skyline.back_substitution(b);
  jsl::map_vector<size_t, double, std::vector>(b, p2d.u, map);
  for (size_t i = 0; i < 5; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
    CHECK(p2d(i, 3) == Approx(p2d.x[i]));
  }

}

TEST_CASE("Case 5 - ADAD, Skyline Incremental 5x3", "[Skyline]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(5, 3);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(M);

  REQUIRE(skyline.heights().size() == 9);
  CHECK(skyline.heights()[0] == 0);
  CHECK(skyline.heights()[1] == 0);
  CHECK(skyline.heights()[2] == 0);
  CHECK(skyline.heights()[3] == 3);
  CHECK(skyline.heights()[4] == 3);
  CHECK(skyline.heights()[5] == 3);
  CHECK(skyline.heights()[6] == 3);
  CHECK(skyline.heights()[7] == 3);
  CHECK(skyline.heights()[8] == 3);

  REQUIRE(skyline.offsets().size() == 9);
  CHECK(skyline.offsets()[0] == 0);
  CHECK(skyline.offsets()[1] == 0);
  CHECK(skyline.offsets()[2] == 0);
  CHECK(skyline.offsets()[3] == 0);
  CHECK(skyline.offsets()[4] == 3);
  CHECK(skyline.offsets()[5] == 6);
  CHECK(skyline.offsets()[6] == 9);
  CHECK(skyline.offsets()[7] == 12);
  CHECK(skyline.offsets()[8] == 15);

  std::vector<double> au{ {-1.0, 0.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0, 
    -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0} };

  REQUIRE(skyline.upper().size() == 18);
  for (size_t i = 0; i < 18; ++i) {
    INFO("The index is " << i);
    CHECK(skyline.upper()[i] == au[i]);
  }

  std::vector<double> v(9);
  jsl::UTDU<size_t, double, std::vector>(9, M, v);

  skyline.utdu();
  for (size_t i = 0; i < 9; i++) {
    INFO("The index is " << i);
    CHECK(skyline.diagonal()[i] == M[i][i]);
  }
  skyline.forward_substitution(b);
  skyline.back_substitution(b);
  jsl::map_vector<size_t, double, std::vector>(b, p2d.u, map);
  for (size_t i = 0; i < 5; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
  }

}

TEST_CASE("Case 5 - ADAD, Skyline", "[Skyline]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(5);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  std::vector<std::vector<double>> M;
  std::vector<double> b;
  std::vector<size_t> map;

  p2d.matrix_system(M, map, b);

  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(M);

  REQUIRE(skyline.heights().size() == 15);
  for (size_t i = 0; i < 3; ++i) {
    CHECK(skyline.heights()[i] == 0);
  }
  for (size_t i = 3; i < 15; ++i) {
    CHECK(skyline.heights()[i] == 3);
  }

  REQUIRE(skyline.offsets().size() == 15);
  for (size_t i = 0; i < 4; ++i) {
    CHECK(skyline.offsets()[i] == 0);
  }
  CHECK(skyline.offsets()[4] == 3);
  CHECK(skyline.offsets()[5] == 6);
  CHECK(skyline.offsets()[6] == 9);
  CHECK(skyline.offsets()[7] == 12);
  CHECK(skyline.offsets()[8] == 15);
  CHECK(skyline.offsets()[9] == 18);
  CHECK(skyline.offsets()[10] == 21);
  CHECK(skyline.offsets()[11] == 24);
  CHECK(skyline.offsets()[12] == 27);
  CHECK(skyline.offsets()[13] == 30);
  CHECK(skyline.offsets()[14] == 33);

  std::vector<double> au{ {-1.0, 0.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, 0.0, 0.0,
    -1.0, 0.0, -1.0, -1.0, 0.0, -1.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0,
    -1.0, 0.0, -1.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0} };

  REQUIRE(skyline.upper().size() == 36);
  for (size_t i = 0; i < 36; ++i) {
    CHECK(skyline.upper()[i] == au[i]);
  }

  skyline.ldlt_solve(b);
  jsl::map_vector<size_t, double, std::vector>(b, p2d.u, map);
  for (size_t i = 0; i < 5; ++i) {
    CHECK(p2d(i, 0) == Approx(p2d.x[i]));
    CHECK(p2d(i, 1) == Approx(p2d.x[i]));
    CHECK(p2d(i, 2) == Approx(p2d.x[i]));
    CHECK(p2d(i, 3) == Approx(p2d.x[i]));
    CHECK(p2d(i, 4) == Approx(p2d.x[i]));
  }
}

TEST_CASE("G&VL Example 4.1.2", "[Skyline]")
{
  std::vector<std::vector<double>> A{ { { 10.0, 20.0, 30.0 }, {20.0, 45.0, 80.0}, {30.0, 80.0, 171.0} } };
  skyline::SymmetricSkyline<size_t, double, std::vector> skyline(A);

  REQUIRE(skyline.heights().size() == 3);
  CHECK(skyline.heights()[0] == 0);
  CHECK(skyline.heights()[1] == 1);
  CHECK(skyline.heights()[2] == 2);

  REQUIRE(skyline.offsets().size() == 3);
  CHECK(skyline.offsets()[0] == 0);
  CHECK(skyline.offsets()[1] == 0);
  CHECK(skyline.offsets()[2] == 1);

  REQUIRE(skyline.minima().size() == 3);
  CHECK(skyline.minima()[0] == 0);
  CHECK(skyline.minima()[0] == 0);
  CHECK(skyline.minima()[0] == 0);

  REQUIRE(skyline.upper().size() == 3);
  CHECK(skyline.upper()[0] == 20.0);
  CHECK(skyline.upper()[1] == 30.0);
  CHECK(skyline.upper()[2] == 80.0);

  REQUIRE(skyline.diagonal().size() == 3);
  CHECK(skyline.diagonal()[0] == 10.0);
  CHECK(skyline.diagonal()[1] == 45.0);
  CHECK(skyline.diagonal()[2] == 171.0);

  skyline.utdu();
  CHECK(skyline.upper()[0] == 2.0);
  CHECK(skyline.upper()[1] == 3.0);
  CHECK(skyline.upper()[2] == 4.0);

  REQUIRE(skyline.diagonal().size() == 3);
  CHECK(skyline.diagonal()[0] == 10.0);
  CHECK(skyline.diagonal()[1] == 5.0);
  CHECK(skyline.diagonal()[2] == 1.0);

  std::vector<double> b{ {0.0, 0.0, 1.0} };

  skyline.forward_substitution(b);
  skyline.back_substitution(b);
  CHECK(b[0] == Approx(5.0));
  CHECK(b[1] == Approx(-4.0));
  CHECK(b[2] == Approx(1.0));
}
