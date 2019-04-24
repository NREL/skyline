// Copyright (c) 2019, Alliance for Sustainable Energy, LLC
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
