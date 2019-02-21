#include "catch.hpp"
#include "../include/poisson2d.hpp"
#include "../dependencies/jsl/jsl.hpp"

#include <iostream>
#include <iomanip>

TEST_CASE("Case 5 - ADAD, GEnxn", "[Poisson2D]")
{
  poisson::Poisson2D<size_t, double, std::vector> p2d(5);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return 1.0; });

  REQUIRE(p2d.ni == 5);
  REQUIRE(p2d.nj == 5);
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

TEST_CASE("Case 1 - DDDD 4x4 f=6xy(1-y)-2x^3, GEnxn", "[Poisson2D]")
{
  size_t N = 4;
  size_t N2 = N * N;
  size_t M2 = (N - 2)*(N - 2);
  poisson::Poisson2D<size_t, double, std::vector> p2d(N);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return y*(1.0 - y); });
  p2d.set_rhs([](double x, double y) { return 6.0*x*y*(1.0 - y) - 2.0*x*x*x; });

  REQUIRE(p2d.ni == N);
  REQUIRE(p2d.nj == N);
  REQUIRE(p2d.x.size() == N2);
  REQUIRE(p2d.y.size() == N2);
  REQUIRE(p2d.u.size() == N2);
  REQUIRE(p2d.f.size() == N2);

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

  REQUIRE(M.size() == M2);
  REQUIRE(map.size() == M2);
  REQUIRE(b.size() == M2);

  std::vector<double> x(M2), z(M2);
  std::vector<size_t> ip(M2);
  jsl::GEnxn<size_t, double, std::vector>(M2, M, x, b, z, ip);

  std::vector<double> spreadsheet{ { 0.0277777777775757, 0.0833333333332323, 0.0277777777776767, 0.0833333333332828} };

  for (size_t i = 0; i < N; ++i) {
    //double x = p2d.x[i];
    //double y = p2d.y[i];
    INFO("Error at index " << i);
    REQUIRE(x[i] == Approx(spreadsheet[i]));
    //CHECK(p2d.u[i] == Approx(y*(1.0 - y)*x*x*x));
  }

  
  jsl::map_vector<size_t, double, std::vector>(x, p2d.u, map);
  for (size_t i = 0; i < N2; ++i) {
    double x = p2d.x[i];
    double y = p2d.y[i];
    //INFO("Error at index " << i);
    //CHECK(p2d.u[i] == Approx(y*(1.0 - y)*x*x*x));
  }

}

TEST_CASE("Case 1 - DDDD 16x16 f=6xy(1-y)-2x^3, GEnxn", "[Poisson2D]")
{
  size_t N = 10;
  size_t N2 = N * N;
  size_t M2 = (N - 2)*(N - 2);
  poisson::Poisson2D<size_t, double, std::vector> p2d(N);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return y * (1.0 - y); });
  p2d.set_rhs([](double x, double y) { return 6.0*x*y*(1.0 - y) - 2.0*x*x*x; });

  REQUIRE(p2d.ni == N);
  REQUIRE(p2d.nj == N);
  REQUIRE(p2d.x.size() == N2);
  REQUIRE(p2d.y.size() == N2);
  REQUIRE(p2d.u.size() == N2);
  REQUIRE(p2d.f.size() == N2);

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

  REQUIRE(M.size() == M2);
  REQUIRE(map.size() == M2);
  REQUIRE(b.size() == M2);

  std::vector<double> x(M2), z(M2);
  std::vector<size_t> ip(M2);
  jsl::GEnxn<size_t, double, std::vector>(M2, M, x, b, z, ip);

  jsl::map_vector<size_t, double, std::vector>(x, p2d.u, map);
  for (size_t i = 0; i < N2; ++i) {
    double x = p2d.x[i];
    double y = p2d.y[i];
    INFO("Error at index " << i);
    //CHECK(p2d.u[i] == Approx(y*(1.0 - y)*x*x*x));
  }

}

TEST_CASE("Case 1 - DDDD 4x4 f=6xy(1-y)-2x^3, GS iteration", "[Poisson2D]")
{
  size_t N = 4;
  size_t N2 = N * N;
  size_t M2 = (N - 2)*(N - 2);
  poisson::Poisson2D<size_t, double, std::vector> p2d(N);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_east([](double y) { return y * (1.0 - y); });
  p2d.set_rhs([](double x, double y) { return 6.0*x*y*(1.0 - y) - 2.0*x*x*x; });

  REQUIRE(p2d.ni == N);
  REQUIRE(p2d.nj == N);
  REQUIRE(p2d.x.size() == N2);
  REQUIRE(p2d.y.size() == N2);
  REQUIRE(p2d.u.size() == N2);
  REQUIRE(p2d.f.size() == N2);

  std::vector<size_t> map;

  auto gs = p2d.gauss_siedel_iterator(map);

  REQUIRE(gs);

  REQUIRE(gs->f.size() == M2);

  //double last = 1.0e6;
  for (int i = 1; i <= 20; ++i) {
    double delta = gs->iterate();
    //REQUIRE(delta <= last);
    //std::cout << i << ' ' << delta << std::endl;
    //last = delta;
    //if (delta == 0.0) {
    //  break;
    //}
  }

  std::vector<double> spreadsheet{ { 0.0277777777775757, 0.0833333333332323, 0.0277777777776767, 0.0833333333332828} };

  jsl::map_vector<size_t, double, std::vector>(gs->u, p2d.u, map);
  for (size_t i = 0; i < N; ++i) {
    //double x = p2d.x[i];
    //double y = p2d.y[i];
    INFO("Error at index " << i);
    REQUIRE(gs->u[i] == Approx(spreadsheet[i]));
    //CHECK(p2d.u[i] == Approx(y*(1.0 - y)*x*x*x));
  }

}


TEST_CASE("Case 1 - DDDD 4x4 f=0, GS iteration", "[Poisson2D]")
{
  size_t N = 4;
  size_t N2 = N * N;
  size_t M2 = (N - 2)*(N - 2);
  poisson::Poisson2D<size_t, double, std::vector> p2d(N);
  p2d.north_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.south_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  p2d.set_south([](double x) { return x; });
  p2d.set_east([](double y) { return 1.0; });
  p2d.set_north([](double x) { return x; });

  REQUIRE(p2d.ni == N);
  REQUIRE(p2d.nj == N);
  REQUIRE(p2d.x.size() == N2);
  REQUIRE(p2d.y.size() == N2);
  REQUIRE(p2d.u.size() == N2);
  REQUIRE(p2d.f.size() == N2);

  std::vector<size_t> map;

  auto gs = p2d.gauss_siedel_iterator(map);

  REQUIRE(gs);

  REQUIRE(gs->f.size() == M2);

  double last = 1.0e6;
  for (int i = 1; i <= 50; ++i) {
    double delta = gs->iterate();
    //REQUIRE(delta <= last);
    //std::cout << i << ' ' << delta << std::endl;
    last = delta;
    //if (delta == 0.0) {
    //  break;
    //}
  }

  jsl::map_vector<size_t, double, std::vector>(gs->u, p2d.u, map);
  for (size_t i = 0; i < N2; ++i) {
    double x = p2d.x[i];
    INFO("Error at index " << i);
    CHECK(p2d.u[i] == Approx(x));
  }

}
