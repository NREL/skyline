#include "catch.hpp"
#include "../src/jsl.hpp"
#include <iostream>

void LUPSolve(std::vector<std::vector<double>> &A, std::vector<double> &b, int N, std::vector<double> &x)
{

  for (int i = 0; i < N; i++) {
    x[i] = b[i];

    for (int k = 0; k < i; k++)
      x[i] -= A[i][k] * x[k];
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++)
      x[i] -= A[i][k] * x[k];

    x[i] = x[i] / A[i][i];
  }
}

TEST_CASE("Small Forward and Backward Substitution, G&VL Example 3.2.2", "[forward_substitution][back_substitution]")
{
  std::vector<std::vector<double>> A{ { { 1.0, 4.0, 7.0 }, {2.0, -3.0, -6.0}, {3.0, 2.0, 1.0} } };
  std::vector<double> b{ {1.0, 1.0, 1.0} };

  jsl::forward_substitution<size_t, double, std::vector>(3, A, b);
  CHECK(b[0] == Approx(1.0));
  CHECK(b[1] == Approx(-1.0));
  CHECK(b[2] == Approx(0.0));

  jsl::back_substitution<size_t, double, std::vector>(3, A, b);

  CHECK(b[0] == Approx(-1.0/3.0));
  CHECK(b[1] == Approx(1.0/3.0));
  CHECK(b[2] == Approx(0.0));

}

TEST_CASE("Small GEnxn", "[GEnxn]")
{
  std::vector<std::vector<double>> A{ { { 10.0, 20.0, 30.0 }, {20.0, 45.0, 80.0}, {30.0, 80.0, 171.0} } };
  std::vector<double> x{ {0.0, 0.0, 0.0} };
  std::vector<double> b{ {0.0, 0.0, 1.0} };
  std::vector<double> z{ {0.0, 0.0, 0.0} };
  std::vector<size_t> ip{ {0, 0, 0} };

  jsl::GEnxn<size_t, double, std::vector>(3, A, x, b, z, ip);

  CHECK(x[0] == Approx(5.0));
  CHECK(x[1] == Approx(-4.0));
  CHECK(x[2] == Approx(1.0));

}

TEST_CASE("Small LDLT, G&VL Example 4.1.2", "[LDLT]")
{
  std::vector<std::vector<double>> A{ { { 10.0, 20.0, 30.0 }, {20.0, 45.0, 80.0}, {30.0, 80.0, 171.0} } };
  std::vector<double> v{ {0.0, 0.0, 0.0} };
  std::vector<double> b{ {0.0, 0.0, 1.0} };
  std::vector<double> D{ {10.0, 5.0, 1.0} };
  
  jsl::LDLT<size_t, double, std::vector>(3, A, v);

  CHECK(A[0][0] == 10.0);
  CHECK(A[0][1] == 20.0);
  CHECK(A[0][2] == 30.0);
  CHECK(A[1][0] == 2.0);
  CHECK(A[1][1] == 5.0);
  CHECK(A[1][2] == 80.0);
  CHECK(A[2][0] == 3.0);
  CHECK(A[2][1] == 4.0);
  CHECK(A[2][2] == 1.0);

  A[0][0] = 1.0;
  A[0][1] = 2.0;
  A[0][2] = 3.0;
  A[1][1] = 1.0;
  A[1][2] = 4.0;
  A[2][2] = 1.0;

  jsl::forward_substitution<size_t, double, std::vector>(3, A, b);

  CHECK(b[0] == Approx(0.0));
  CHECK(b[1] == Approx(0.0));
  CHECK(b[2] == Approx(1.0));

  // Account for the diagonal
  for (size_t i = 0; i < 3; ++i) {
    b[i] /= D[i];
  }

  jsl::back_substitution<size_t, double, std::vector>(3, A, b);

  CHECK(b[0] == Approx(5.0));
  CHECK(b[1] == Approx(-4.0));
  CHECK(b[2] == Approx(1.0));

}

TEST_CASE("Small UTDU, G&VL Example 4.1.2", "[UTDU]")
{
  std::vector<std::vector<double>> A{ { { 10.0, 20.0, 30.0 }, {20.0, 45.0, 80.0}, {30.0, 80.0, 171.0} } };
  std::vector<double> v{ {0.0, 0.0, 0.0} };

  jsl::UTDU<size_t, double, std::vector>(3, A, v);

  CHECK(A[0][0] == 10.0);
  CHECK(A[0][1] == 2.0);
  CHECK(A[0][2] == 3.0);
  CHECK(A[1][0] == 20.0);
  CHECK(A[1][1] == 5.0);
  CHECK(A[1][2] == 4.0);
  CHECK(A[2][0] == 30.0);
  CHECK(A[2][1] == 80.0);
  CHECK(A[2][2] == 1.0);

}
