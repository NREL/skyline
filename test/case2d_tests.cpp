#include "catch.hpp"
#include "../include/poisson2d.hpp"

TEST_CASE( "Case 1 - DDDD", "[Case2D]" ) {
  auto case1 = poisson::Case2D<int>::diagnose(4, 3, poisson::BoundaryCondition::Dirichlet, poisson::BoundaryCondition::Dirichlet,
    poisson::BoundaryCondition::Dirichlet, poisson::BoundaryCondition::Dirichlet);
  REQUIRE(case1);
  CHECK(case1->type == poisson::Type::DDDD);
  CHECK(case1->ni == 4);
  CHECK(case1->nj == 3);
  CHECK(case1->mi == 2);
  CHECK(case1->mj == 1);
  CHECK(case1->stride == 2);
  CHECK(case1->start == 5);
  // Test other stuff
}

TEST_CASE("Case 5 - DADA", "[Case2D]")
{
  auto case5 = poisson::Case2D<int>::diagnose(5, 8, poisson::BoundaryCondition::Dirichlet, poisson::BoundaryCondition::Adiabatic,
    poisson::BoundaryCondition::Dirichlet, poisson::BoundaryCondition::Adiabatic);
  REQUIRE(case5);
  CHECK(case5->type == poisson::Type::DADA);
  CHECK(case5->rotation == poisson::Rotation::CCW0);
  CHECK(case5->ni == 5);
  CHECK(case5->nj == 8);
  CHECK(case5->mi == 5);
  CHECK(case5->mj == 6);
  CHECK(case5->stride == 0);
  CHECK(case5->start == 5);
  // Test other stuff
}

TEST_CASE("Case 5 - ADAD", "[Case2D]")
{
  auto case5 = poisson::Case2D<int>::diagnose(5, 8, poisson::BoundaryCondition::Adiabatic, poisson::BoundaryCondition::Dirichlet,
    poisson::BoundaryCondition::Adiabatic, poisson::BoundaryCondition::Dirichlet);
  REQUIRE(case5);
  CHECK(case5->type == poisson::Type::DADA);
  CHECK(case5->rotation == poisson::Rotation::CCW90);
  CHECK(case5->ni == 5);
  CHECK(case5->nj == 8);
  CHECK(case5->mi == 3);
  CHECK(case5->mj == 8);
  CHECK(case5->stride == 2);
  CHECK(case5->start == 1);
  // Test other stuff
}
