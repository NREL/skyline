#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <functional>
#include <iomanip>
#include <optional>
#include "skyline.hpp"
#include "poisson2d.hpp"
#include "jsl.hpp"

int main()
{
  poisson::Poisson2D<size_t, double, std::vector> box0(7);
  box0.north_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  box0.south_boundary_condition = poisson::BoundaryCondition::Adiabatic;
  box0.east_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  box0.west_boundary_condition = poisson::BoundaryCondition::Dirichlet;
  box0.set_west([](double y) { return 1.0; });

  //int iter_max = 500;
  //for (int i = 0; i < iter_max; ++i) {
  //  double delta = box0.gauss_seidel_iteration();
    //box.adiabatic_north();
    //box.adiabatic_south();
  //  std::cout << i + 1 << ' ' << delta << std::endl;
  //  if (delta < 1.0e-17) {
  //    break;
  //  }
  //}

  for (size_t i = 0; i < 7; ++i) {
    std::cout << i << ' ' << box0(i, 0) << ' ' << box0(i, 3) << std::endl;
  }
  std::cout << std::endl;

  //UnitSquareGrid<size_t, double, std::vector> box1(7);
  //box1.north_boundary_condition = BoundaryCondition::Adiabatic;
  //box1.south_boundary_condition = BoundaryCondition::Adiabatic;
  //box1.east_boundary_condition = BoundaryCondition::Dirichlet;
  //box1.west_boundary_condition = BoundaryCondition::Dirichlet;
  //box1.set_west([](double y) { return 1.0; });

  std::vector<std::vector<double>> A;
  std::vector<double> f;
  std::vector<size_t> x_map;

  size_t N = box0.matrix_system(A, x_map, f);

  for (auto v : x_map) {
    std::cout << v << std::endl;
  }

  for (size_t i = 0; i < A.size(); ++i) {
    for (size_t j = 0; j < A.size(); ++j) {
      std::cout << std::setw(2) << A[i][j] << ' ';
    }
    std::cout << "| " << f[i] << std::endl;
  }


  std::cout << jsl::is_symmetric<size_t, double, std::vector>(A) << std::endl;

  // Solve the system with Gaussian elimination
  std::vector<size_t> ip(N);
  std::vector<double> xge(N);
  std::vector<double> temp(N);
  for (size_t i = 0; i < N; ++i) {
    ip[i] = 0;
    temp[i] = 0.0;
    xge[i] = 0.0;
  }

  jsl::GEnxn<size_t, double, std::vector>(N, A, xge, f, temp, ip);
  jsl::map_vector<size_t, double, std::vector>(xge, box0.u, x_map);

  std::cout << std::endl;
  for (size_t i = 0; i < 7; ++i) {
    std::cout << i << ' ' << box0(i, 0) << ' ' << box0(i, 3) << std::endl;
  }
  std::cout << std::endl;

  /*
  int i = 0;
  for (auto &v : skyline.heights()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;

  i = 0;
  for (auto &v : skyline.offsets()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;

  i = 0;
  for (auto &v : skyline.upper()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;
  */
  exit(EXIT_SUCCESS);

  /*
  std::vector<std::vector<double>> M{ {{4.0, 1.0, 1.0, 0.0, 0.0},
                                       {1.0, 6.0, 0.0, 1.0, 0.0},
                                       {1.0, 0.0, 5.0, 0.0, 0.0},
                                       {0.0, 1.0, 0.0, 7.0, 1.0},
                                       {0.0, 0.0, 0.0, 1.0, 9.0}} };

  //skyline::IndexSolver<size_t, double, std::vector> skyline(M);

  i = 0;
  for (auto &v : skyline.heights()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;

  i = 0;
  for (auto &v : skyline.offsets()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;

  i = 0;
  for (auto &v : skyline.upper()) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }
  std::cout << std::endl;

  for (auto &a : M) {
    for (auto &v : a) {
      std::cout << v << ' ';
    }
    std::cout << std::endl;
  }

  std::vector<int> ik;
  ik.push_back(0);
  for (auto &v : skyline.offsets()) {
    ik.push_back((int)v + 1);
  }
  std::cout << std::endl;
  i = 0;
  for (auto &v : ik) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }

  std::vector<double> au, ad, x;
  au.push_back(0.0);
  ad.push_back(0.0);
  x.push_back(0.0);
  for (auto v : skyline.upper()) {
    au.push_back(v);
  }
  for (auto v : skyline.diagonal()) {
    ad.push_back(v);
    x.push_back(1.0);
  }

  skyline::FACSKYmod(au, ad, au, ik, 5, 0);

  i = 0;
  std::cout << std::endl;
  for (auto &v : ad) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }

  i = 0;
  std::cout << std::endl;
  for (auto &v : au) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }

  skyline::SLVSKYmod(au, ad, au, x, ik, 5, 0);

  i = 0;
  std::cout << std::endl;
  for (auto &v : x) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }

  // Check
  std::vector<double> b{ {0.0, 0.0, 0.0, 0.0, 0.0} };
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < 5; ++j) {
      b[i] += M[i][j] * x[j + 1];
    }
  }

  i = 0;
  std::cout << std::endl;
  for (auto &v : b) {
    std::cout << i << ' ' << v << std::endl;
    ++i;
  }

  return EXIT_SUCCESS;
  */
}
