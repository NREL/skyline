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
#ifndef POISSON2D_HPP
#define POISSON2D_HPP

#include <optional>
#include <iostream>

namespace poisson {

enum class BoundaryCondition { Dirichlet, Adiabatic, Periodic};
enum class Type { DDDD=1, DDDA, DDAA, DAAA, DADA, AAAA, DPDP, APAP, DPAP, PPPP};
enum class Rotation { CCW0 = 0, CCW90 = 90, CCW180 = 180, CCW270 = 270 };

template<typename I, typename R, template <typename ...> typename V> struct GaussSiedelIterator
{

  GaussSiedelIterator(I ni, I nj, V<R> &f, BoundaryCondition N = BoundaryCondition::Dirichlet,
    BoundaryCondition E = BoundaryCondition::Dirichlet, BoundaryCondition S = BoundaryCondition::Dirichlet,
    BoundaryCondition W = BoundaryCondition::Dirichlet) : ni(ni), nj(nj), north_boundary_condition(N), east_boundary_condition(E),
    south_boundary_condition(S), west_boundary_condition(W), f(f)
  {
    u.resize(f.size());
    for (I i = 0; i < f.size(); ++i) {
      u[i] = 0;
    }
  }

  R iterate()
  {
    R delta = 0.0;
    I ij = 0;
    // First row
    if (west_boundary_condition == BoundaryCondition::Adiabatic) {
      u[0] = u[1];
      delta = std::max(delta, std::abs(u[0] - u[1]));
      ij = 1;
    } else if (west_boundary_condition == BoundaryCondition::Dirichlet 
      && south_boundary_condition == BoundaryCondition::Dirichlet) {
      R last = u[0];
      u[0] = 0.25*(f[0] + u[1] + u[ni]);
      delta = std::abs(last - u[0]);
      ij = 1;
    }

    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (; ij < ni-1; ++ij) {
        delta = std::max(delta, std::abs(u[ij] - u[ij + ni]));
        u[ij] = u[ij + ni];
      }
    } else if (south_boundary_condition == BoundaryCondition::Dirichlet) {
      for (; ij < ni - 1; ++ij) {
        R last = u[ij];
        u[ij] = 0.25*(f[ij] + u[ij + 1] + u[ij - 1] + u[ij + ni]);
        delta = std::max(delta, std::abs(last - u[ij]));
      }
    }

    if (east_boundary_condition == BoundaryCondition::Adiabatic) {
      u[ij] = u[ij-1];
      delta = std::max(delta, std::abs(u[ij] - u[ij-1]));
      ++ij;
    } else if (east_boundary_condition == BoundaryCondition::Dirichlet
      && south_boundary_condition == BoundaryCondition::Dirichlet) {
      R last = u[ij];
      u[ij] = 0.25*(f[ij] + u[ij-1] + u[ij+ni]);
      delta = std::abs(last - u[ij]);
      ++ij;
    }

    // Main body
    for (I j = 1; j < nj-1; ++j) {
      // First
      if (west_boundary_condition == BoundaryCondition::Adiabatic) {
      } else {
        R last = u[ij];
        u[ij] = 0.25*(f[ij] + u[ij + 1] + u[ij + ni] + u[ij - ni]);
        delta = std::max(delta, std::abs(last - u[ij]));
      }
      ++ij;
      for (I i = 1; i < ni-1; ++i) {
        R last = u[ij];
        u[ij] = 0.25*(f[ij] + u[ij - 1] + u[ij + 1] + u[ij + ni] + u[ij - ni]);
        delta = std::max(delta, std::abs(last - u[ij]));
        ++ij;
      }
      // Last
      if (east_boundary_condition == BoundaryCondition::Adiabatic) {
      } else {
        R last = u[ij];
        u[ij] = 0.25*(f[ij] + u[ij - 1] + u[ij + ni] + u[ij - ni]);
        delta = std::max(delta, std::abs(last - u[ij]));
      }
      ++ij;
    }
    
    // Last row
    if (west_boundary_condition == BoundaryCondition::Adiabatic) {
      delta = std::max(delta, std::abs(u[ij] - u[ij + 1]));
      u[ij] = u[ij+1];
    } else if (west_boundary_condition == BoundaryCondition::Dirichlet
      && south_boundary_condition == BoundaryCondition::Dirichlet) {
      R last = u[ij];
      u[ij] = 0.25*(f[0] + u[ij+1] + u[ij-ni]);
      delta = std::max(delta, std::abs(last - u[ij]));
    }
    ++ij;

    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (; ij < ni - 1; ++ij) {
        delta = std::max(delta, std::abs(u[ij] - u[ij - ni]));
        u[ij] = u[ij - ni];
      }
    } else if (south_boundary_condition == BoundaryCondition::Dirichlet) {
      for (; ij < ni - 1; ++ij) {
        R last = u[ij];
        u[ij] = 0.25*(f[ij] + u[ij + 1] + u[ij - 1] + u[ij + ni]);
        delta = std::max(delta, std::abs(last - u[ij]));
      }
    }

    if (east_boundary_condition == BoundaryCondition::Adiabatic) {
      u[ij] = u[ij - 1];
      delta = std::max(delta, std::abs(u[ij] - u[ij - 1]));
    } else if (east_boundary_condition == BoundaryCondition::Dirichlet
      && north_boundary_condition == BoundaryCondition::Dirichlet) {
      R last = u[ij];
      u[ij] = 0.25*(f[ij] + u[ij - 1] + u[ij - ni]);
      delta = std::abs(last - u[ij]);
    }
    return delta;
  }

  I ni, nj;
  BoundaryCondition north_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition east_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition south_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition west_boundary_condition = BoundaryCondition::Dirichlet;
  V<R> u, f;
};

template <typename I> struct Case2D
{
  Case2D(Type type, Rotation rotation, I ni, I nj, I mi, I mj, I start, I stride) : type(type), rotation(rotation), ni(ni), nj(nj), mi(mi), mj(mj),
    start(start), stride(stride)
  {}

  static std::optional<Case2D> diagnose(I ni, I nj, BoundaryCondition N, BoundaryCondition E, BoundaryCondition S, BoundaryCondition W)
  {
    int nD{ 0 };
    int nP{ 0 };
    if (N == BoundaryCondition::Dirichlet) {
      ++nD;
    } else if(N == BoundaryCondition::Periodic) {
      ++nP;
    }
    if (E == BoundaryCondition::Dirichlet) {
      ++nD;
    } else if(E == BoundaryCondition::Periodic) {
      ++nP;
    }
    if (S == BoundaryCondition::Dirichlet) {
      ++nD;
    } else if(S == BoundaryCondition::Periodic) {
      ++nP;
    }
    if (W == BoundaryCondition::Dirichlet) {
      ++nD;
    } else if(W == BoundaryCondition::Periodic) {
      ++nP;
    }
    if (nP > 0) {
      // To do at some point
    } else {
      // Non-periodic cases
      if (nD == 4) {
        // DDDD
        return Case2D(Type::DDDD, Rotation::CCW0, ni, nj, ni - 2, nj - 2, ni + 1, 2);
      } else if(nD == 3) {
        // DDDA
        Rotation rotation{ Rotation::CCW0 };
        I mi{ ni - 1 };
        I mj{ nj - 2 };
        I start{ ni };
        I stride{ 1 };
        if (N == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW90;
          mi = ni - 2;
          mj = nj - 1;
          start = 1;
          stride = 2;
        } else if(E == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW180;
          mi = ni - 1;
          mj = nj - 2;
          start = ni + 1;
          stride = 1;
        } else if(S == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW270;
          mi = ni - 2;
          mj = nj - 1;
          start = ni + 1;
          stride = 2;
        }
        return Case2D(Type::DDDA, rotation, ni, nj, mi, mj, start, stride);
      } else if (nD == 1) {
        // DAAA
        Rotation rotation{ Rotation::CCW0 };
        I mi{ ni };
        I mj{ nj - 1 };
        I start{ 0 };
        I stride{ 0 };
        if (W == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW90;
          mi = ni - 1;
          mj = nj;
          start = 1;
          stride = 1;
        } else if (S == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW180;
          mi = ni;
          mj = nj - 1;
          start = ni;
          stride = 1;
        } else if (E == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW270;
          mi = ni - 1;
          mj = nj;
          start = 0;
          stride = 1;
        }
        return Case2D(Type::DAAA, rotation, ni, nj, mi, mj, start, stride);
      } else if (nD == 0) {
        // AAAA, To do at some point
      } else {
        // DDAA or DADA
        if (N == BoundaryCondition::Dirichlet) {
          // DDAA(0), DAAD(90), DADA(0)
          if (S == BoundaryCondition::Dirichlet) {
            // DADA
            return Case2D(Type::DADA, Rotation::CCW0, ni, nj, ni, nj-2, ni, 0);
          } else if (W == BoundaryCondition::Dirichlet) {
            // DAAD
            return Case2D(Type::DDAA, Rotation::CCW90, ni, nj, ni-1, nj-1, 1, 1);
          }
          // DDAA
          return Case2D(Type::DDAA, Rotation::CCW0, ni, nj, ni - 1, nj - 1, 0, 1);
        } else {
          // AADD(180), ADDA(270), ADAD(90)
          if (E == BoundaryCondition::Adiabatic) {
            // AADD
            return Case2D(Type::DDAA, Rotation::CCW180, ni, nj, ni - 1, nj - 1, ni + 1, 1);
          } else if (S == BoundaryCondition::Dirichlet) {
            // ADDA
            return Case2D(Type::DDAA, Rotation::CCW270, ni, nj, ni - 1, nj - 1, ni, 1);
          }
          // ADAD
          return Case2D(Type::DADA, Rotation::CCW90, ni, nj, ni - 2, nj, 1, 2);
        }
      }
    }
    return {};
  }

  Type type{ Type::DDDD };
  Rotation rotation{ Rotation::CCW0 };
  I ni{ 0 }, nj{ 0 }, mi{ 0 }, mj{ 0 };
  I start{ 0 }, stride{ 0 };
};

template <typename I, typename R, template <typename ...> typename V> struct Poisson2D
{

  Poisson2D(I ni, I nj) : ni(std::max(ni, (I)3)), nj(std::max(nj, (I)3))
  {
    I N2 = ni * nj;
    x.resize(N2);
    y.resize(N2);
    u.resize(N2);
    f.resize(N2);
    for (I i = 0; i < N2; ++i) {
      u[i] = 0.0;
      f[i] = 0.0;
    }
    R deltax = 1.0 / (R)(ni - 1);
    R deltay = 1.0 / (R)(nj - 1);
    I ij = 0;
    for (I j = 0; j < nj; ++j) {
      R yy = deltay * j;
      for (I i = 0; i < ni; ++i) {
        y[ij] = yy;
        x[ij] = deltax * i;
        ++ij;
      }
    }
  }

  Poisson2D(I n) : Poisson2D(n,n)
  {}

  R operator()(I i, I j)
  {
    return u[i + j * ni];
  }

  void set_west(std::function<R(R)> fcn)
  {
    I ij = 0;
    for (I j = 0; j < nj; ++j) {
      u[ij] = fcn(y[ij]);
      ij += ni;
    }
  }

  void set_east(std::function<R(R)> fcn)
  {
    I ij = ni-1;
    for (I j = 0; j < nj; ++j) {
      u[ij] = fcn(y[ij]);
      ij += ni;
    }
  }

  void set_south(std::function<R(R)> fcn)
  {
    I ij = 0;
    for (I j = 0; j < nj; ++j) {
      u[ij] = fcn(x[ij]);
      ++ij;
    }
  }

  void set_north(std::function<R(R)> fcn)
  {
    I ij = nj*(ni - 1);
    for (I j = 0; j < nj; ++j) {
      u[ij] = fcn(x[ij]);
      ++ij;
    }
  }

  void set_rhs(std::function<R(R,R)> fcn)
  {
    I ij = nj * (ni - 1);
    for (I j = 0; j < nj; ++j) {
      f[ij] = fcn(x[ij],y[ij]);
      ++ij;
    }
  }

  const I ni, nj;
  V<R> x, y, u, f;
  BoundaryCondition north_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition east_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition south_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition west_boundary_condition = BoundaryCondition::Dirichlet;

  std::optional<GaussSiedelIterator<I,R,V>> gauss_siedel_iterator(V<I> &x_map, std::ostream *out = nullptr)
  {
    auto twod = Case2D<I>::diagnose(ni, nj, north_boundary_condition,
      east_boundary_condition, south_boundary_condition, west_boundary_condition);
    if (!twod) {
      return {};
    }

    if (out != nullptr) {
      *out << " Type: " << (int)(twod->type) << std::endl;
      *out << "  Rot: " << (int)(twod->rotation) << std::endl;
      *out << "   ni: " << twod->ni << std::endl;
      *out << "   nj: " << twod->nj << std::endl;
      *out << "   mi: " << twod->mi << std::endl;
      *out << "   mj: " << twod->mj << std::endl;
      *out << "Start: " << twod->start << std::endl;
      *out << "Shift: " << twod->stride << std::endl;
    }

    I N = twod->mi*twod->mj;
    V<R> newf(N);
    x_map.resize(N);
    I ij = twod->start;
    I mij = 0;
    for (I j = 0; j < twod->mj; ++j) {
      for (I i = 0; i < twod->mi; ++i) {
        x_map[mij] = ij;
        newf[mij] = f[ij];
        ++mij;
        ++ij;
      }
      ij += twod->stride;
    }

    // Modify for boundary conditions
    // Handle W border
    ij = twod->start;
    mij = 0;
    if (west_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I j = 0; j < twod->mj; ++j) {
        newf[mij] += u[ij - 1];
        mij += twod->mi;
        ij += ni;
      }
    } else {
      mij += twod->mj*twod->mi;
      ij += twod->mj*ni;
    }
    // Handle N border
    mij = (twod->mj - 1)*twod->mi;
    ij = x_map[mij];
    if (north_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I i = 0; i < twod->mi; ++i) {
        newf[mij] += u[ij + ni];
        ++mij;
        ++ij;
      }
    }

    // Handle S border
    ij = twod->start;
    mij = 0;
    if (south_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I i = 0; i < twod->mi; ++i) {
        newf[mij] += u[ij - ni];
        ++mij;
        ++ij;
      }
    }
    // Handle E border
    mij = twod->mi-1;
    ij = x_map[mij];
    if (east_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I j = 0; j < twod->mj; ++j) {
        newf[mij] += u[ij + 1];
        mij += twod->mi;
        ij += ni;
      }
    }
    
    return GaussSiedelIterator<I, R, V>(twod->mi, twod->mj, newf, north_boundary_condition, east_boundary_condition,
      south_boundary_condition, west_boundary_condition);
    
  }

  I matrix_system(V<V<R>> &M, V<I> &x_map, V<R> &b, std::ostream *out = nullptr)
  {
    auto twod = Case2D<I>::diagnose(ni, nj, north_boundary_condition,
      east_boundary_condition, south_boundary_condition, west_boundary_condition);
    if (!twod) {
      return 0;
    }

    if (out != nullptr) {
      *out << " Type: " << (int)(twod->type) << std::endl;
      *out << "  Rot: " << (int)(twod->rotation) << std::endl;
      *out << "   ni: " << twod->ni << std::endl;
      *out << "   nj: " << twod->nj << std::endl;
      *out << "   mi: " << twod->mi << std::endl;
      *out << "   mj: " << twod->mj << std::endl;
      *out << "Start: " << twod->start << std::endl;
      *out << "Shift: " << twod->stride << std::endl;
    }

    I N = twod->mi*twod->mj;

    x_map.resize(N);
    I ij = twod->start;
    I mij = 0;
    for (I j = 0; j < twod->mj; ++j) {
      for (I i = 0; i < twod->mi; ++i) {
        x_map[mij] = ij;
        ++mij;
        ++ij;
      }
      ij += twod->stride;
    }

    // Initialize the matrix equations
    b.resize(twod->mi*twod->mj);
    M.resize(twod->mi*twod->mj);
    for (I i = 0; i < N; ++i) {
      M[i].resize(twod->mi*twod->mj);
      b[i] = f[x_map[i]];
      for (I j = 0; j < N; ++j) {
        M[i][j] = 0.0;
      }
    }

    // Set up the equations
    ij = twod->start;
    mij = 0;
    // j = 0
    // Handle SW corner
    if (west_boundary_condition == BoundaryCondition::Dirichlet &&
      south_boundary_condition == BoundaryCondition::Dirichlet) {
      b[mij] += u[ij - 1] + u[ij - ni];
      M[mij][mij] = 4.0;
      M[mij][mij + twod->mi] = -1.0;
      M[mij][mij + 1] = -1.0;
      ++ij;
      ++mij;
    }
    // Now the rest
    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = 0; i < twod->mi; ++i) {
        M[mij][mij] = 1.0;
        M[mij][mij + twod->mi] = -1.0;
        ++ij;
        ++mij;
      }
    } else if (south_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I i = 1; i < twod->mi-1; ++i) {
        b[mij] += u[ij - ni];
        M[mij][mij] = 4.0;
        M[mij][mij + twod->mi] = -1.0;
        M[mij][mij - 1] = -1.0;
        M[mij][mij + 1] = -1.0;
        ++ij;
        ++mij;
      }
    }
    // Handle the SE corner
    if (east_boundary_condition == BoundaryCondition::Dirichlet &&
      south_boundary_condition == BoundaryCondition::Dirichlet) {
      b[mij] += u[ij + 1] + u[ij - ni];
      M[mij][mij] = 4.0;
      M[mij][mij + twod->mi] = -1.0;
      M[mij][mij - 1] = -1.0;
      ++ij;
      ++mij;
    }
    ij += twod->stride;
    // 1 <= j < mj-1
    for (I j = 1; j < twod->mj-1; ++j) {
      // i = 0
      if (west_boundary_condition == BoundaryCondition::Dirichlet) {
        b[mij] += u[ij - 1];
        M[mij][mij] = 4.0;
        M[mij][mij - twod->mi] = -1.0;
        M[mij][mij + twod->mi] = -1.0;
        M[mij][mij + 1] = -1.0;
      } else {
        // others go here
      }
      ++ij;
      ++mij;
      for (I i = 1; i < twod->mi-1; ++i) {
        M[mij][mij] = 4.0;
        M[mij][mij - twod->mi] = -1.0;
        M[mij][mij + twod->mi] = -1.0;
        M[mij][mij - 1] = -1.0;
        M[mij][mij + 1] = -1.0;
        ++mij;
        ++ij;
      }
      // i = mi-1
      if (west_boundary_condition == BoundaryCondition::Dirichlet) {
        b[mij] += u[ij + 1];
        M[mij][mij] = 4.0;
        M[mij][mij - twod->mi] = -1.0;
        M[mij][mij + twod->mi] = -1.0;
        M[mij][mij - 1] = -1.0;
      } else {
        // others go here
      }
      ++mij;
      ++ij;
      ij += twod->stride;
    }
    // j = mj-1
    // Handle the NW corner
    if (west_boundary_condition == BoundaryCondition::Dirichlet &&
      north_boundary_condition == BoundaryCondition::Dirichlet) {
      b[mij] += u[ij - 1] + u[ij + ni];
      M[mij][mij] = 4.0;
      M[mij][mij - twod->mi] = -1.0;
      M[mij][mij + 1] = -1.0;
      ++ij;
      ++mij;
    }
    // Now the rest
    if (north_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = 0; i < twod->mi; ++i) {
        M[mij][mij] = 1.0;
        M[mij][mij - twod->mi] = -1.0;
        ++ij;
        ++mij;
      }
    } else if (north_boundary_condition == BoundaryCondition::Dirichlet) {
      for (I i = 1; i < twod->mi - 1; ++i) {
        b[mij] += u[ij + ni];
        M[mij][mij] = 4.0;
        M[mij][mij - twod->mi] = -1.0;
        M[mij][mij - 1] = -1.0;
        M[mij][mij + 1] = -1.0;
        ++ij;
        ++mij;
      }
    }
    // Handle the NE corner
    if (east_boundary_condition == BoundaryCondition::Dirichlet &&
      north_boundary_condition == BoundaryCondition::Dirichlet) {
      b[mij] += u[ij + 1] + u[ij + ni];
      M[mij][mij] = 4.0;
      M[mij][mij - twod->mi] = -1.0;
      M[mij][mij - 1] = -1.0;
      ++ij;
      ++mij;
    }
    return N;
  }

};

}

#endif
