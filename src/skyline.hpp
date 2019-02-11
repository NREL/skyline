#ifndef SKYLINE_HPP
#define SKYLINE_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <optional>

template <typename T> class OneArray
{
public:
  OneArray(int n, T default = 0)
  {
    for (int i = 0; i < n; ++i) {
      m_v.push_back(default);
    }
  }

  T operator[](int i) const
  {
    return m_v[i - 1];
  }
  T operator()(int i) const
  {
    return m_v[i - 1];
  }
  T& operator[](int i)
  {
    return m_v[i - 1];
  }
  T& operator()(int i)
  {
    return m_v[i - 1];
  }

private:
  std::vector<T> m_v;
};


namespace skyline {

template <typename I, typename R, template <typename ...> typename V> class SymmetricSkyline
{
public:

  SymmetricSkyline(V<V<R>> &M)
  {
    I n = M.size();
    for (auto &v : M) {
      n = std::min(n, v.size());
    }
    m_ih.resize(n);
    for (I i = 0; i < n; i++) {
      m_ih[i] = i;
      for (I j = 0; j < i; j++) {
        if (M[i][j] != 0.0) {
          break;
        }
        --m_ih[i];
      }
    }
    auto sum = std::accumulate(m_ih.begin(), m_ih.end(), (I)0);

    m_ik.resize(n);
    m_im.resize(n);
    // Convert heights to column offsets.
    m_ik[0] = 0;
    m_im[0] = 0;
    for (I k = 1; k < n; ++k) {
      m_ik[k] = m_ik[k - 1] + m_ih[k - 1];
      m_im[k] = k - m_ih[k];
    }

    // Copy into the au/d array
    m_ad.resize(n);
    m_au.resize(sum);
    I count = 0;
    for (I i = 0; i < n; i++) {
      m_ad[i] = M[i][i];
      I j;
      for (j = 0; j < i; j++) {
        if (M[i][j] != 0.0) {
          break;
        }
      }
      for (; j < i; j++) {
        m_au[count] = M[i][j];
        ++count;
      }
    }

    m_v.resize(n);
    m_n = n;
  }

  V<I> offsets() const
  {
    return m_ik;
  }

  V<I> heights() const
  {
    return m_ih;
  }

  V<I> minima() const
  {
    return m_im;
  }

  V<R> diagonal() const
  {
    return m_ad;
  }

  V<R> upper() const
  {
    return m_au;
  }

  V<R> lower() const
  {
    return m_au;
  }

  std::optional<I> index(I i, I j) const
  {
    if (m_im[j] <= i) {
      return m_ik[j] + i - m_im[j];
    }
    return {};
  }

  void utdu()
  {
    // j = 0, nothing much to do
    for (I k = 1; k < m_n; ++k) {
      if (m_im[k] == 0) {
        I ij = m_ik[k];
        m_au[ij] /= m_ad[0];
      }
    }
    // Now for the rest
    for (I j = 1; j < m_n; ++j) {
      // Compute v
      for (I i = 0; i < m_im[j]; ++i) {
        m_v[i] = 0.0;
      }
      for (I i = m_im[j]; i < j; ++i) {
        m_v[i] = m_au[m_ik[j] + i - m_im[j]] * m_ad[i];
      }
      // Compute the diagonal term
      R value = 0.0;
      for (I i = m_im[j]; i < j; ++i) {
        value += m_au[m_ik[j] + i - m_im[j]] * m_v[i];
      }
      m_ad[j] -= value;
      // Compute the rest of the row
      for (I k = j + 1; k < m_n; ++k) {
        if (m_im[k] <= j) {
          value = 0.0;
          for (I i = m_im[k]; i < j; ++i) {
            I ij = m_ik[k] + i - m_im[k];
            value += m_au[ij] * m_v[i];
          }
          I ij = m_ik[j] + k - m_im[j];
          m_au[ij] = (m_au[ij] - value) / m_ad[j];
        }
      }
    }
  }

  void forward_substitution(V<R> &b) const
  {
    // Solve Lz=b (Dy=z, Ux=y)
    for (I i = 1; i < m_n; ++i) {
      R value = 0.0;
      for (I k = m_im[i]; k < i; ++k) {
        I ij = m_ik[i] + k - m_im[i];
        value += m_au[ij] * b[k];
      }
      b[i] -= value;
    }
  }

  void back_substitution(V<R> &z) const
  {
    // Account for the diagonal first (invert Dy=z)
    for (I j = 0; j < m_n; ++j) {
      z[j] /= m_ad[j];
    }
    // Solve Ux=y
    for (I j = m_n - 1; j > 0; --j) {
      for (I k = m_im[j]; k < j; ++k) {
        I ij = m_ik[j] + k - m_im[j];
        z[k] -= z[j] * m_au[ij];
      }
    }
  }

  void ldlt_solve(V<R> &b)
  {
    utdu();
    forward_substitution(b);
    back_substitution(b);
  }

private:

  void factorize()
  {

    // SUBROUTINE INFORMATION:
    //       AUTHOR         George Walton
    //       DATE WRITTEN   Extracted from AIRNET
    //       MODIFIED       Lixing Gu, 2/1/04
    //                      Revised the subroutine to meet E+ needs
    //       MODIFIED       Lixing Gu, 6/8/05
    //       RE-ENGINEERED  This subroutine is revised from FACSKY developed by George Walton, NIST

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine performs L-U factorization of a skyline ordered matrix, [A]
    // The algorithm has been restructured for clarity.
    // Note dependence on compiler for optimizing the inner do loops.

    // METHODOLOGY EMPLOYED:
    //     L-U factorization of a skyline ordered matrix, [A], used for
    //     solution of simultaneous linear algebraic equations [A] * X = B.
    //     No pivoting!  No scaling!  No warnings!!!
    //     Related routines:  SLVSKY, SETSKY, FILSKY.

    // REFERENCES:
    //     Algorithm is described in "The Finite Element Method Displayed",
    //     by G. Dhatt and G. Touzot, John Wiley & Sons, New York, 1984.

    // USE STATEMENTS:
    // na

    // Argument array dimensioning
    //AU.dim(IK(NetworkNumOfNodes + 1));
    //AD.dim(NetworkNumOfNodes);
    //AL.dim(IK(NetworkNumOfNodes + 1) - 1);
    //IK.dim(NetworkNumOfNodes + 1);

    // Locals
    // SUBROUTINE ARGUMENT DEFINITIONS:
    // noel, GNU says the AU is indexed above its upper bound
    // REAL(r64), INTENT(INOUT) :: AU(IK(NetworkNumOfNodes+1)-1) ! the upper triangle of [A] before and after factoring

    // SUBROUTINE PARAMETER DEFINITIONS:
    // na

    // INTERFACE BLOCK SPECIFICATIONS
    // na

    // DERIVED TYPE DEFINITIONS
    // na

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    //Real64 T1;
    //Real64 T2;
    //Real64 SDOT;

    I NEQ{ m_ik.size() };

    // FLOW:
    m_ad[0] = 1.0 / m_ad[0];
    I JHK = 1;
    for (I k = 1; k < NEQ; ++k) {
      Real64 SUMD = 0.0;
      I JHK1 = m_ik[k + 1];
      I LHK = JHK1 - JHK;
      if (LHK > 0) {
        I LHK1 = LHK - 1;
        I IMIN = k - LHK1;
        I IMIN1 = IMIN - 1;
        if (!symmetric) {
          m_al[JHK] *= m_ad[IMIN1];
        }
        if (LHK1 != 0) {
          I JHJ = m_ik[IMIN];
          if (symmetric) {
            for (I j = 0; j < LHK1; ++j) {
              I JHJ1 = m_ik[IMIN + j];
              I IC = std::min(j, JHJ1 - JHJ);
              if (IC > 0) {
                real_t SDOT = 0.0;
                for (I i = 0; i <= IC - 1; ++i) {
                  SDOT += m_au[JHJ1 - IC + i] * m_au[JHK + j - IC + i];
                }
                m_au[JHK + j] -= SDOT;
              }
              JHJ = JHJ1;
            }
          } else {
            /*
            for (j = 1; j <= LHK1; ++j) {
              JHJ1 = IK[IMIN + j];
              IC = std::min(j, JHJ1 - JHJ);
              SDOT = 0.0;
              if (IC > 0) {
                for (i = 0; i <= IC - 1; ++i) {
                  SDOT += AL[JHJ1 - IC + i] * AU[JHK + j - IC + i];
                }
                AU[JHK + j] -= SDOT;
                SDOT = 0.0;
                for (i = 0; i <= IC - 1; ++i) {
                  SDOT += AU[JHJ1 - IC + i] * AL[JHK + j - IC + i];
                }
              }
              AL[JHK + j] = (AL[JHK + j] - SDOT) * AD[IMIN1 + j];
              JHJ = JHJ1;
            }
            */
          }
        }
        if (symmetric) {
          for (I i = 0; i <= LHK1; ++i) {
            real_t T1 = m_au[JHK + i];
            real_t T2 = T1 * m_ad[IMIN1 + i];
            m_au[JHK + i] = T2;
            SUMD += T1 * T2;
          }
        } else {
          /*
          for (i = 0; i <= LHK1; ++i) {
            SUMD += AU[JHK + i] * AL[JHK + i];
          }
          */
        }
      }
      if (m_ad[k] == SUMD) {
        std::cerr << "The denominator used in L-U factorization is equal to 0.0 at index " << k << std::endl;
        exit(EXIT_FAILURE);
        //ShowSevereError("AirflowNetworkSolver: L-U factorization in Subroutine FACSKY.");
        //ShowContinueError("The denominator used in L-U factorizationis equal to 0.0 at node = " + AirflowNetworkNodeData(k).Name + '.');
        //ShowContinueError(
        //  "One possible cause is that this node may not be connected directly, or indirectly via airflow network connections ");
        //ShowContinueError(
        //  "(e.g., AirflowNetwork:Multizone:SurfaceCrack, AirflowNetwork:Multizone:Component:SimpleOpening, etc.), to an external");
        //ShowContinueError("node (AirflowNetwork:MultiZone:Surface).");
        //ShowContinueError("Please send your input file and weather file to EnergyPlus support/development team for further investigation.");
        //ShowFatalError("Preceding condition causes termination.");
      }
      m_ad[k] = 1.0 / (m_ad[k] - SUMD);
      JHK = JHK1;
    }
  }

  I m_n;
  V<I> m_ik; // Index offsets to top of skylines
  V<I> m_ih; // Height of each skyline
  V<I> m_im; // Minimum row, or top of skyline
  V<R> m_au; // Upper triangular part of matrix
  V<R> m_ad; // Diagonal of matrix
  V<R> m_v;  // Temporary used in solution
};


template <typename I, typename R, template <typename ...> typename V> class IndexSolver
{
public:

  IndexSolver(V<V<R>> &M) : symmetric(true)
  {
    I n = M.size();
    for (auto &v : M) {
      n = std::min(n, v.size());
    }
    m_ih.resize(n);
    for (I i = 0; i < n; i++) {
      m_ih[i] = i;
      for (I j = 0; j < i; j++) {
        if (M[i][j] != 0.0) {
          break;
        }
        --m_ih[i];
      }
    }
    auto sum = std::accumulate(m_ih.begin(), m_ih.end(), (I)0);
    setsky();
    // Copy into the au/d array
    m_ad.resize(n);
    m_au.resize(sum);
    I count = 0;
    for (I i = 0; i < n; i++) {
      m_ad[i] = M[i][i];
      I j;
      for (j = 0; j < i; j++) {
        if (M[i][j] != 0.0) {
          break;
        }
      }
      for (; j < i; j++) {
        m_au[count] = M[i][j];
        ++count;
      }
    }
  }

  V<I> offsets() const
  {
    return m_ik;
  }

  V<I> heights() const
  {
    return m_ih;
  }

  V<R> diagonal() const
  {
    return m_ad;
  }

  V<R> upper() const
  {
    return m_au;
  }

  V<R> lower() const
  {
    return m_au;
  }

  const bool symmetric;

private:

  void setsky()
  {
    // SUBROUTINE INFORMATION:
    //       AUTHOR         George Walton
    //       DATE WRITTEN   1998
    //       MODIFIED       Feb. 2006 (L. Gu) to meet requirements of AirflowNetwork
    //       RE-ENGINEERED  na

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine sets up the "IK" array describing the sparse matrix [A] in skyline
    //     form by using the location matrix.

    // METHODOLOGY EMPLOYED:
    // na

    // REFERENCES:
    // AIRNET

    // USE STATEMENTS:
    // na

    // Locals
    // SUBROUTINE ARGUMENT DEFINITIONS:
    // na

    // SUBROUTINE PARAMETER DEFINITIONS:
    // na

    // INTERFACE BLOCK SPECIFICATIONS
    // na

    // DERIVED TYPE DEFINITIONS
    // na

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    // IK(K) - pointer to the top of column/row "K".

    // FLOW:
    // Initialize "IK".
    //for (i = 1; i <= NetworkNumOfNodes + 1; ++i) {
    //  IK[i] = 0;
    //}
    // Determine column heights.
    //for (M = 1; M <= NetworkNumOfLinks; ++M) {
    //  j = AirflowNetworkLinkageData(M).NodeNums[1];
    //  if (j == 0) continue;
    //  L = ID(j);
    //  i = AirflowNetworkLinkageData(M).NodeNums[0];
    //  k = ID(i);
    //  N1 = std::abs(L - k);
    //  N2 = max(k, L);
    //  IK(N2) = max(IK(N2), N1);
    //}
    // Set things up
    I N = m_ih.size();
    //m_ih = h;
    //m_au.resize(N);
    //m_ad.resize(N);
    m_ik.resize(N);
    // Convert heights to column addresses.
    m_ik[0] = 0;
    for (I k = 1; k < N; ++k) {
      m_ik[k] = m_ik[k - 1] + m_ih[k - 1];
    }
  }


  void factorize()
  {

    // SUBROUTINE INFORMATION:
    //       AUTHOR         George Walton
    //       DATE WRITTEN   Extracted from AIRNET
    //       MODIFIED       Lixing Gu, 2/1/04
    //                      Revised the subroutine to meet E+ needs
    //       MODIFIED       Lixing Gu, 6/8/05
    //       RE-ENGINEERED  This subroutine is revised from FACSKY developed by George Walton, NIST

    // PURPOSE OF THIS SUBROUTINE:
    // This subroutine performs L-U factorization of a skyline ordered matrix, [A]
    // The algorithm has been restructured for clarity.
    // Note dependence on compiler for optimizing the inner do loops.

    // METHODOLOGY EMPLOYED:
    //     L-U factorization of a skyline ordered matrix, [A], used for
    //     solution of simultaneous linear algebraic equations [A] * X = B.
    //     No pivoting!  No scaling!  No warnings!!!
    //     Related routines:  SLVSKY, SETSKY, FILSKY.

    // REFERENCES:
    //     Algorithm is described in "The Finite Element Method Displayed",
    //     by G. Dhatt and G. Touzot, John Wiley & Sons, New York, 1984.

    // USE STATEMENTS:
    // na

    // Argument array dimensioning
    //AU.dim(IK(NetworkNumOfNodes + 1));
    //AD.dim(NetworkNumOfNodes);
    //AL.dim(IK(NetworkNumOfNodes + 1) - 1);
    //IK.dim(NetworkNumOfNodes + 1);

    // Locals
    // SUBROUTINE ARGUMENT DEFINITIONS:
    // noel, GNU says the AU is indexed above its upper bound
    // REAL(r64), INTENT(INOUT) :: AU(IK(NetworkNumOfNodes+1)-1) ! the upper triangle of [A] before and after factoring

    // SUBROUTINE PARAMETER DEFINITIONS:
    // na

    // INTERFACE BLOCK SPECIFICATIONS
    // na

    // DERIVED TYPE DEFINITIONS
    // na

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    //Real64 T1;
    //Real64 T2;
    //Real64 SDOT;

    I NEQ{ m_ik.size() };

    // FLOW:
    m_ad[0] = 1.0 / m_ad[0];
    I JHK = 1;
    for (I k = 1; k < NEQ; ++k) {
      Real64 SUMD = 0.0;
      I JHK1 = m_ik[k + 1];
      I LHK = JHK1 - JHK;
      if (LHK > 0) {
        I LHK1 = LHK - 1;
        I IMIN = k - LHK1;
        I IMIN1 = IMIN - 1;
        if (!symmetric) {
          m_al[JHK] *= m_ad[IMIN1];
        }
        if (LHK1 != 0) {
          I JHJ = m_ik[IMIN];
          if (symmetric) {
            for (I j = 0; j < LHK1; ++j) {
              I JHJ1 = m_ik[IMIN + j];
              I IC = std::min(j, JHJ1 - JHJ);
              if (IC > 0) {
                real_t SDOT = 0.0;
                for (I i = 0; i <= IC - 1; ++i) {
                  SDOT += m_au[JHJ1 - IC + i] * m_au[JHK + j - IC + i];
                }
                m_au[JHK + j] -= SDOT;
              }
              JHJ = JHJ1;
            }
          } else {
            /*
            for (j = 1; j <= LHK1; ++j) {
              JHJ1 = IK[IMIN + j];
              IC = std::min(j, JHJ1 - JHJ);
              SDOT = 0.0;
              if (IC > 0) {
                for (i = 0; i <= IC - 1; ++i) {
                  SDOT += AL[JHJ1 - IC + i] * AU[JHK + j - IC + i];
                }
                AU[JHK + j] -= SDOT;
                SDOT = 0.0;
                for (i = 0; i <= IC - 1; ++i) {
                  SDOT += AU[JHJ1 - IC + i] * AL[JHK + j - IC + i];
                }
              }
              AL[JHK + j] = (AL[JHK + j] - SDOT) * AD[IMIN1 + j];
              JHJ = JHJ1;
            }
            */
          }
        }
        if (symmetric) {
          for (I i = 0; i <= LHK1; ++i) {
            real_t T1 = m_au[JHK + i];
            real_t T2 = T1 * m_ad[IMIN1 + i];
            m_au[JHK + i] = T2;
            SUMD += T1 * T2;
          }
        } else {
          /*
          for (i = 0; i <= LHK1; ++i) {
            SUMD += AU[JHK + i] * AL[JHK + i];
          }
          */
        }
      }
      if (m_ad[k] == SUMD) {
        std::cerr << "The denominator used in L-U factorization is equal to 0.0 at index " << k << std::endl;
        exit(EXIT_FAILURE);
        //ShowSevereError("AirflowNetworkSolver: L-U factorization in Subroutine FACSKY.");
        //ShowContinueError("The denominator used in L-U factorizationis equal to 0.0 at node = " + AirflowNetworkNodeData(k).Name + '.');
        //ShowContinueError(
        //  "One possible cause is that this node may not be connected directly, or indirectly via airflow network connections ");
        //ShowContinueError(
        //  "(e.g., AirflowNetwork:Multizone:SurfaceCrack, AirflowNetwork:Multizone:Component:SimpleOpening, etc.), to an external");
        //ShowContinueError("node (AirflowNetwork:MultiZone:Surface).");
        //ShowContinueError("Please send your input file and weather file to EnergyPlus support/development team for further investigation.");
        //ShowFatalError("Preceding condition causes termination.");
      }
      m_ad[k] = 1.0 / (m_ad[k] - SUMD);
      JHK = JHK1;
    }
  }


  V<I> m_ik;
  V<I> m_ih;
  V<R> m_au;
  V<R> m_ad;
  V<R> m_al;
};

void SETSKY(std::vector<int> &IK);

void FACSKYmod(std::vector<double> &AU,   // the upper triangle of [A] before and after factoring
  std::vector<double> &AD,   // the main diagonal of [A] before and after factoring
  std::vector<double> &AL,   // the lower triangle of [A] before and after factoring
  std::vector<int> const &IK, // pointer to the top of column/row "K"
  int const NEQ,        // number of equations
  int const NSYM        // symmetry:  0 = symmetric matrix, 1 = non-symmetric
);

void SLVSKYmod(std::vector<double> const &AU, // the upper triangle of [A] before and after factoring
  std::vector<double> const &AD, // the main diagonal of [A] before and after factoring
  std::vector<double> const &AL, // the lower triangle of [A] before and after factoring
  std::vector<double> &B,        // "B" vector (input); "X" vector (output).
  std::vector<int> const &IK,     // pointer to the top of column/row "K"
  int const NEQ,            // number of equations
  int const NSYM            // symmetry:  0 = symmetric matrix, 1 = non-symmetric
);


void FACSKY(OneArray<double> AU,   // the upper triangle of [A] before and after factoring
  OneArray<double> AD,   // the main diagonal of [A] before and after factoring
  OneArray<double> AL,   // the lower triangle of [A] before and after factoring
  OneArray<int> const IK, // pointer to the top of column/row "K"
  int const NEQ,        // number of equations
  int const NSYM        // symmetry:  0 = symmetric matrix, 1 = non-symmetric
);

void SLVSKY(OneArray<double> const AU, // the upper triangle of [A] before and after factoring
  OneArray<double> const AD, // the main diagonal of [A] before and after factoring
  OneArray<double> const AL, // the lower triangle of [A] before and after factoring
  OneArray<double> B,        // "B" vector (input); "X" vector (output).
  OneArray<int> const IK,     // pointer to the top of column/row "K"
  int const NEQ,            // number of equations
  int const NSYM            // symmetry:  0 = symmetric matrix, 1 = non-symmetric
);

}

#endif // !SKYLINE_HPP
