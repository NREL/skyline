#ifndef POISSON2D_HPP
#define POISSON2D_HPP

#include <optional>
#include <iostream>

namespace poisson {

enum class BoundaryCondition { Dirichlet, Adiabatic, Periodic};
enum class Type { DDDD=1, DDDA, DDAA, DAAA, DADA, AAAA, DPDP, APAP, DPAP, PPPP};
enum class Rotation { CCW0 = 0, CCW90 = 90, CCW180 = 180, CCW270 = 270 };

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

  //Poisson2D(I ni, I nj)
  //{
  //}

  Poisson2D(I n) : N(std::max(n, (I)4))
  {
    I N2 = N * N;
    x.resize(N2);
    y.resize(N2);
    u.resize(N2);
    for (I i = 0; i < N2; ++i) {
      u[i] = 0.0;
    }
    R delta = 1.0 / (R)(N-1);
    I ij = 0;
    for (I j = 0; j < N; ++j) {
      R yy = delta * j;
      for (I i = 0; i < N; ++i) {
        y[ij] = yy;
        x[ij] = delta * i;
        ++ij;
      }
    }
  }

  R operator()(I i, I j)
  {
    return u[i + j * N];
  }

  void set_west(std::function<R(R)> f)
  {
    I ij = 0;
    for (I j = 0; j < N; ++j) {
      u[ij] = f(y[ij]);
      ij += N;
    }
  }

  void set_east(std::function<R(R)> f)
  {
    I ij = N-1;
    for (I j = 0; j < N; ++j) {
      u[ij] = f(y[ij]);
      ij += N;
    }
  }

  void set_south(std::function<R(R)> f)
  {
    I ij = 0;
    for (I j = 0; j < N; ++j) {
      u[ij] = f(x[ij]);
      ++ij;
    }
  }

  void set_north(std::function<R(R)> f)
  {
    I ij = N*(N - 1);
    for (I j = 0; j < N; ++j) {
      u[ij] = f(x[ij]);
      ++ij;
    }
  }

  const I N;
  V<R> x, y, u;
  BoundaryCondition north_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition south_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition east_boundary_condition = BoundaryCondition::Dirichlet;
  BoundaryCondition west_boundary_condition = BoundaryCondition::Dirichlet;


  R gauss_seidel_iteration()
  {
    R delta = 0.0;
    I ij = N + 1;
    I i1 = 1;
    I i2 = N - 1;
    I ijshift = 2;

    I j1 = 1;
    I j2 = N-1;
    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = i1; i < i2; ++i) {
        delta = std::max(delta, std::abs(u[i] - u[i + N]));
        u[i] = u[i + N];
      }
    }
    for (I j = j1; j < j2; ++j) {
      for (I i = i1; i < i2; ++i) {
        R last = u[ij];
        u[ij] = 0.25*(u[ij - 1] + u[ij + 1] + u[ij + N] + u[ij - N]);
        delta = std::max(delta, std::abs(last - u[ij]));
        ++ij;
      }
      ij += ijshift;
    }
    if (north_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = N*N-N; i < N*N-1; ++i) {
        delta = std::max(delta, std::abs(u[i] - u[i - N]));
        u[i] = u[i - N];
      }
    }
    return delta;
  }

  void adiabatic_south()
  {
    I ij = 0;
    for (I j = 0; j < N; ++j) {
      u[ij] = u[ij+N];
      ++ij;
    }
  }

  void adiabatic_north()
  {
    I ij = N * (N - 1);
    for (I j = 0; j < N; ++j) {
      u[ij] = u[ij-N];
      ++ij;
    }
  }

  I matrix_system(V<V<R>> &M, V<I> &x_map, V<R> &b, std::ostream *out = nullptr)
  {
    auto twod = Case2D<I>::diagnose(N, N, north_boundary_condition,
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
      b[i] = 0.0;
      for (I j = 0; j < N; ++j) {
        M[i][j] = 0.0;
      }
    }

    // Set up the equations
    ij = twod->start;
    mij = 0;
    // j = 0
    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = 0; i < twod->mi; ++i) {
        M[mij][mij] = 1.0;
        M[mij][mij + twod->mi] = -1.0;
        ++ij;
        ++mij;
      }
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
    if (north_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = 0; i < twod->mi; ++i) {
        M[mij][mij] = 1.0;
        M[mij][mij - twod->mi] = -1.0;
        ++ij;
        ++mij;
      }
    }

    return N;
  }

};

}

#endif
