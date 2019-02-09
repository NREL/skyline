#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <functional>
#include <iomanip>
#include <optional>
#include "skyline.hpp"

typedef double Real64;

enum class BoundaryCondition { Dirichlet, Adiabatic, Periodic};
enum class TypeN { Type1 = 1, Type2, Type3, Type4, Type5, Type6, Type7, Type8, Type9, Type10 };
enum class Type { DDDD=1, DDDA, DDAA, DAAA, DADA, AAAA, DPDP, APAP, DPAP, PPPP};
enum class Rotation { CCW0 = 0, CCW90 = 90, CCW180 = 180, CCW270 = 270 };

template <typename I> struct TwoDimensionalCase
{
  TwoDimensionalCase(Type type, Rotation rotation, I ni, I nj, I mi, I mj, I nstart, I nshift) : type(type), rotation(rotation), ni(ni), nj(nj), mi(mi), mj(mj),
    nstart(nstart), nshift(nshift)
  {}

  static std::optional<TwoDimensionalCase> diagnose(I ni, I nj, BoundaryCondition N, BoundaryCondition E, BoundaryCondition S, BoundaryCondition W)
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
        return TwoDimensionalCase(Type::DDDD, Rotation::CCW0, ni, nj, ni - 1, nj - 1, ni + 1, 2);
      } else if(nD == 3) {
        // DDDA
        Rotation rotation{ Rotation::CCW0 };
        I mi{ ni - 1 };
        I mj{ nj - 2 };
        I nstart{ ni };
        I nshift{ 1 };
        if (N == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW90;
          mi = ni - 2;
          mj = nj - 1;
          nstart = 1;
          nshift = 2;
        } else if(E == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW180;
          mi = ni - 1;
          mj = nj - 2;
          nstart = ni + 1;
          nshift = 1;
        } else if(S == BoundaryCondition::Adiabatic) {
          rotation = Rotation::CCW270;
          mi = ni - 2;
          mj = nj - 1;
          nstart = ni + 1;
          nshift = 2;
        }
        return TwoDimensionalCase(Type::DDDA, rotation, ni, nj, mi, mj, nstart, nshift);
      } else if (nD == 1) {
        // DAAA
        Rotation rotation{ Rotation::CCW0 };
        I mi{ ni };
        I mj{ nj - 1 };
        I nstart{ 0 };
        I nshift{ 0 };
        if (W == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW90;
          mi = ni - 1;
          mj = nj;
          nstart = 1;
          nshift = 1;
        } else if (S == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW180;
          mi = ni;
          mj = nj - 1;
          nstart = ni;
          nshift = 1;
        } else if (E == BoundaryCondition::Dirichlet) {
          rotation = Rotation::CCW270;
          mi = ni - 1;
          mj = nj;
          nstart = 0;
          nshift = 1;
        }
        return TwoDimensionalCase(Type::DAAA, rotation, ni, nj, mi, mj, nstart, nshift);
      } else if (nD == 0) {
        // AAAA, To do at some point
      } else {
        // DDAA or DADA
        if (N == BoundaryCondition::Dirichlet) {
          // DDAA(0), DAAD(90), DADA(0)
          if (S == BoundaryCondition::Dirichlet) {
            // DADA
            return TwoDimensionalCase(Type::DADA, Rotation::CCW0, ni, nj, ni, nj-2, ni, 0);
          } else if (W == BoundaryCondition::Dirichlet) {
            // DAAD
            return TwoDimensionalCase(Type::DDAA, Rotation::CCW90, ni, nj, ni-1, nj-1, 1, 1);
          }
          // DDAA
          return TwoDimensionalCase(Type::DDAA, Rotation::CCW0, ni, nj, ni - 1, nj - 1, 0, 1);
        } else {
          // AADD(180), ADDA(270), ADAD(90)
          if (E == BoundaryCondition::Adiabatic) {
            // AADD
            return TwoDimensionalCase(Type::DDAA, Rotation::CCW180, ni, nj, ni - 1, nj - 1, ni + 1, 1);
          } else if (S == BoundaryCondition::Dirichlet) {
            // ADDA
            return TwoDimensionalCase(Type::DDAA, Rotation::CCW270, ni, nj, ni - 1, nj - 1, ni, 1);
          }
          // ADAD
          return TwoDimensionalCase(Type::DADA, Rotation::CCW90, ni, nj, ni - 2, nj, 1, 2);
        }
      }
    }
    return {};
  }

  Type type{ Type::Type1 };
  Rotation rotation{ Rotation::CCW0 };
  I ni{ 0 }, nj{ 0 }, mi{ 0 }, mj{ 0 };
  I nstart{ 0 }, nshift{ 0 };
};

template <typename I, typename R, template <typename ...> typename V> bool is_symmetric(V<V<R>> &M)
{
  I N = M.size();
  for (I j = 0; j < N; ++j) {
    for (I i = j+1; i < N; ++i) {
      if (M[i][j] != M[j][i]) {
        std::cout << i << ' ' << j << ' ' << M[i][j] << " != " << M[j][i] << std::endl;
        return false;
      }
    }
  }
  return true;
}

template <typename I, typename R, template <typename ...> typename V> void map_vector(V<R> &from, V<R> &to, V<I> map)
{
  I N = from.size();
  for (I i = 0; i < N; ++i) {
    to[map[i]] = from[i];
  }
}

template <typename I, typename R, template <typename ...> typename V> struct UnitSquareGrid
{

  UnitSquareGrid(I n) : N(std::max(n, (I)4))
  {
    I N2 = N * N;
    x.resize(N2);
    y.resize(N2);
    u.resize(N2);
    for (I i = 0; i < N2; ++i) {
      u[i] = 0.0;
    }
    R delta = 1.0 / (R)N;
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

  void set_east(std::function<R(R)> f)
  {
    I ij = 0;
    for (I j = 0; j < N; ++j) {
      u[ij] = f(y[ij]);
      ij += N;
    }
  }

  void set_west(std::function<R(R)> f)
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

  I matrix_system(V<V<R>> &M, V<I> &x_map, V<R> &b)
  {
    auto twod = TwoDimensionalCase<I>::diagnose(N, N, north_boundary_condition,
      east_boundary_condition, south_boundary_condition, west_boundary_condition);
    if (!twod) {
      return 0;
    }

    std::cout << " Type: " << (int)(twod->type) << std::endl;
    std::cout << "  Rot: " << (int)(twod->rotation) << std::endl;
    std::cout << "   ni: " << twod->ni << std::endl;
    std::cout << "   nj: " << twod->nj << std::endl;
    std::cout << "   mi: " << twod->mi << std::endl;
    std::cout << "   mj: " << twod->mj << std::endl;
    std::cout << "Start: " << twod->nstart << std::endl;
    std::cout << "Shift: " << twod->nshift << std::endl;
    //TwoDimensionalCase<I> twod{opttwod.get()};

    I N = twod->mi*twod->mj;

    x_map.resize(N);
    I ij = twod->nstart;
    I mij = 0;
    for (I j = 0; j < twod->mj; ++j) {
      for (I i = 0; i < twod->mi; ++i) {
        x_map[mij] = ij;
        ++mij;
        ++ij;
      }
      ij += twod->nshift;
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
    ij = twod->nstart;
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
    ij += twod->nshift;
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
      ij += twod->nshift;
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
    /*
    I ni = N;
    I nj = N;
    I ij = 0;
    I i1 = 0;
    I i2 = N;
    I ijshift = 0;
    I j1 = 0;
    I j2 = N;
    if (north_boundary_condition == BoundaryCondition::Dirichlet) {
      --nj;
      --j2;
    }
    if (south_boundary_condition == BoundaryCondition::Dirichlet) {
      --nj;
      ++j1;
      ij += N;
    }
    if (east_boundary_condition == BoundaryCondition::Dirichlet) {
      --ni;
      ++ij;
      ++i1;
      ++ijshift;
    }
    if (west_boundary_condition == BoundaryCondition::Dirichlet) {
      --ni;
      --i2;
      ++ijshift;
    }

    I ijstart = ij;

    std::cout << "System is " << ni << " by " << nj << std::endl;
    std::cout << "i1 : " << i1 << std::endl;
    std::cout << "i2 : " << i2 << std::endl;
    std::cout << "j1 : " << j1 << std::endl;
    std::cout << "j2 : " << j2 << std::endl;
    std::cout << "ijstart : " << ijstart << std::endl;

    b.resize(ni*nj);
    for (I j = 0; j < ni*nj; ++j) {
      b[j] = 0.0;
    }

    x_map.resize(ni*nj);
    I ijmap = 0;
    for (I j = j1; j < j2; ++j) {
      for (I i = i1; i < i2; ++i) {
        x_map[ijmap] = ij;
        ++ijmap;
        ++ij;
      }
      ij += ijshift;
    }

    b.resize(ni*nj);
    M.resize(ni*nj);
    for (I i = 0; i < ni*nj; ++i) {
      M[i].resize(ni*nj);
      for (I j = 0; j < ni*nj; ++j) {
        M[i][j] = 0.0;
      }
    }

    // Handle Boundary conditions
    ijmap = 0;
    if (south_boundary_condition == BoundaryCondition::Adiabatic) {

    }

    ijstart = ni;
    ijshift = 2;

    // North
    ij = ni * nj - ni;
    if (north_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = i1; i < i2; ++i) {
        M[ij][ij] = 1.0;
        M[ij][ij - ni] = -1.0;
        ++ij;
      }
      --j2;
    } //else

    // South
    ij = 0;
    if (south_boundary_condition == BoundaryCondition::Adiabatic) {
      for (I i = i1; i < i2; ++i) {
        M[ij][ij] = 1.0;
        M[ij][ij + ni] = -1.0;
        ++ij;
      }
      ++j1;
    } // else
    // West
    if (west_boundary_condition == BoundaryCondition::Adiabatic) {
      // TODO
    } else {
      ijmap = j1 * ni;
      for (I j = j1; j < j2; j++) {
        ij = x_map[ijmap];
        b[ijmap] += u[ij-1];
        M[ijmap][ijmap] = 4.0;
        M[ijmap][ijmap - ni] = -1.0;
        M[ijmap][ijmap + ni] = -1.0;
        M[ijmap][ijmap + 1] = -1.0;
      }
      ++i1;
      ++ijstart;
    }
    // East
    if (east_boundary_condition == BoundaryCondition::Adiabatic) {
      // TODO
    } else {
      ijmap = (j1+1) * ni - 1;
      for (I j = j1; j < j2; j++) {
        ij = x_map[ijmap];
        b[ijmap] += u[ij+1];
        M[ijmap][ijmap] = 4.0;
        M[ijmap][ijmap - ni] = -1.0;
        M[ijmap][ijmap + ni] = -1.0;
        M[ijmap][ijmap - 1] = -1.0;
      }
      --i2;
    }


    std::cout << "i1 : " << i1 << std::endl;
    std::cout << "i2 : " << i2 << std::endl;
    std::cout << "j1 : " << j1 << std::endl;
    std::cout << "j2 : " << j2 << std::endl;
    std::cout << "ijstart : " << ijstart << std::endl;

    ij = ijstart;
    for (I j = j1; j < j2; ++j) {
      for (I i = i1; i < i2; ++i) {
        M[ij][ij] = 4.0;
        M[ij][ij - ni] = -1.0;
        M[ij][ij + ni] = -1.0;
        M[ij][ij - 1] = -1.0;
        M[ij][ij + 1] = -1.0;
        ++ij;
      }
      ij += ijshift;
    }
    */
    return true;
  }

};

/*****************************************************************************/
/*                                                                           */
/*   GEnxn - Do Gaussian elimination on a n by n system of equations.        */
/*                                                                           */
/*   Arguments:                                                              */
/*      int_t n ------------ number of rows and columns.                     */
/*      real_t A[][] ------- LHS matrix.                                     */
/*      real_t x[] --------- solution vector.                                */
/*      real_t b[] --------- RHS vector.                                     */
/*      real_t *z ---------- real work vector at least n in length.          */
/*      int_t *ip ---------- integer work vector at least n in length.       */
/*                                                                           */
/*   This routine does Gaussian Elimination on a nxn matrix with partial     */
/*   pivoting.  It probably should be replaced with something better.  A is  */
/*   modified and b is not.  It is a descendant of the F77 version, but is   */
/*   way better because templating.                                          */
/*                                                                           */
/*****************************************************************************/
template <typename I, typename R, template <typename ...> typename V> void GEnxn(I n, V<V<R>> &A, V<R> &x, V<R> &b,
  V<R> &z, V<I> &ip)
{
  R lik, big;
  I i, j, k, itemp;
  /*
    Initialize the setup and copy the RHS
  */
  for (j = 0; j < n; j++) {
    ip[j] = j;
    z[j] = b[j];
  }
  /*
    Forward elimination
  */
  for (k = 0; k < n - 1; k++) {
    big = std::abs(A[ip[k]][k]);
    for (i = k + 1; i < n; i++) {
      if (std::abs(A[ip[i]][k]) > big) {
        itemp = ip[i];
        ip[i] = ip[k];
        ip[k] = itemp;
      }
    }
    /* write(*,*)'Using row: ',ip(k) */
    for (i = k + 1; i < n; i++) {
      lik = A[ip[i]][k] / A[ip[k]][k];
      for (j = k + 1; j < n; j++) {
        A[ip[i]][j] = A[ip[i]][j] - lik * A[ip[k]][j];
      }
      z[ip[i]] = z[ip[i]] - lik * z[ip[k]];
    }
  }

  /*      i=0
c      call print3x3(A,i)
c      write(*,*)(z(k),k=1,3)
c      write(*,*)(ip(j),j=1,3)
  */
  /*
    Back substitute, with tricky looping in case I is unsigned
  */
  k = n;
  do {
    --k;
    for (i = k + 1; i < n; i++) {
      z[ip[k]] = z[ip[k]] - A[ip[k]][i] * z[ip[i]];
    }
    z[ip[k]] = z[ip[k]] / A[ip[k]][k];
  } while (k > 0);
  /*
    Reorder the solution as we copy it to the output
  */
  for (i = 0; i < n; i++) {
    x[i] = z[ip[i]];
  }
  return;
}

int main()
{
  UnitSquareGrid<size_t, double, std::vector> box0(7);
  box0.north_boundary_condition = BoundaryCondition::Adiabatic;
  box0.south_boundary_condition = BoundaryCondition::Adiabatic;
  box0.east_boundary_condition = BoundaryCondition::Dirichlet;
  box0.west_boundary_condition = BoundaryCondition::Dirichlet;
  box0.set_west([](double y) { return 1.0; });

  int iter_max = 500;
  for (int i = 0; i < iter_max; ++i) {
    double delta = box0.gauss_seidel_iteration();
    //box.adiabatic_north();
    //box.adiabatic_south();
    std::cout << i + 1 << ' ' << delta << std::endl;
    if (delta < 1.0e-17) {
      break;
    }
  }

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

  skyline::IndexSolver<size_t, double, std::vector> skyline(A);

  std::cout << is_symmetric<size_t, double, std::vector>(A) << std::endl;

  // Solve the system with Gaussian elimination
  std::vector<size_t> ip(N);
  std::vector<double> xge(N);
  std::vector<double> temp(N);
  for (size_t i = 0; i < N; ++i) {
    ip[i] = 0;
    temp[i] = 0.0;
    xge[i] = 0.0;
  }

  GEnxn<size_t, double, std::vector>(N, A, xge, f, temp, ip);
  map_vector<size_t, double, std::vector>(xge, box0.u, x_map);

  std::cout << std::endl;
  for (size_t i = 0; i < 7; ++i) {
    std::cout << i << ' ' << box0(i, 0) << ' ' << box0(i, 3) << std::endl;
  }
  std::cout << std::endl;


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

  exit(EXIT_SUCCESS);


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
  std::vector<Real64> b{ {0.0, 0.0, 0.0, 0.0, 0.0} };
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
}
