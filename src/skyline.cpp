#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "skyline.hpp"

typedef double Real64;

namespace skyline {

void SETSKY(std::vector<int> &IK)
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
  int i;
  int j;
  int k;
  //int L;
  //int M;
  //int N1;
  //int N2;

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
  // Convert heights to column addresses.
  j = IK[1];
  IK[1] = 1;
  for (k = 1; k <= IK.size(); ++k) {
    i = IK[k + 1];
    IK[k + 1] = IK[k] + j;
    j = i;
  }
}

void FACSKYmod(std::vector<Real64> &AU,   // the upper triangle of [A] before and after factoring
  std::vector<Real64> &AD,   // the main diagonal of [A] before and after factoring
  std::vector<Real64> &AL,   // the lower triangle of [A] before and after factoring
  std::vector<int> const &IK, // pointer to the top of column/row "K"
  int const NEQ,        // number of equations
  int const NSYM        // symmetry:  0 = symmetric matrix, 1 = non-symmetric
)
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
  int JHK;
  int JHK1;
  int LHK;
  int LHK1;
  int IMIN;
  int IMIN1;
  int JHJ;
  int JHJ1;
  int IC;
  int i;
  int j;
  int k;
  Real64 T1;
  Real64 T2;
  Real64 SDOT;
  Real64 SUMD;

  // FLOW:
  AD[1] = 1.0 / AD[1];
  JHK = 1;
  for (k = 2; k <= NEQ; ++k) {
    SUMD = 0.0;
    JHK1 = IK[k];
    LHK = JHK1 - JHK;
    if (LHK > 0) {
      LHK1 = LHK - 1;
      IMIN = k - LHK1;
      IMIN1 = IMIN - 1;
      if (NSYM == 1) AL[JHK] *= AD[IMIN1];
      if (LHK1 != 0) {
        JHJ = IK[IMIN];
        if (NSYM == 0) {
          for (j = 1; j <= LHK1; ++j) {
            JHJ1 = IK[IMIN + j];
            IC = std::min(j, JHJ1 - JHJ);
            if (IC > 0) {
              SDOT = 0.0;
              for (i = 0; i <= IC - 1; ++i) {
                SDOT += AU[JHJ1 - IC + i] * AU[JHK + j - IC + i];
              }
              AU[JHK + j] -= SDOT;
            }
            JHJ = JHJ1;
          }
        } else {
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
        }
      }
      if (NSYM == 0) {
        for (i = 0; i <= LHK1; ++i) {
          T1 = AU[JHK + i];
          T2 = T1 * AD[IMIN1 + i];
          AU[JHK + i] = T2;
          SUMD += T1 * T2;
        }
      } else {
        for (i = 0; i <= LHK1; ++i) {
          SUMD += AU[JHK + i] * AL[JHK + i];
        }
      }
    }
    if (AD[k] - SUMD == 0.0) {
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
    AD[k] = 1.0 / (AD[k] - SUMD);
    JHK = JHK1;
  }
}

void SLVSKYmod(std::vector<Real64> const &AU, // the upper triangle of [A] before and after factoring
  std::vector<Real64> const &AD, // the main diagonal of [A] before and after factoring
  std::vector<Real64> const &AL, // the lower triangle of [A] before and after factoring
  std::vector<Real64> &B,        // "B" vector (input); "X" vector (output).
  std::vector<int> const &IK,     // pointer to the top of column/row "K"
  int const NEQ,            // number of equations
  int const NSYM            // symmetry:  0 = symmetric matrix, 1 = non-symmetric
)
{

  // SUBROUTINE INFORMATION:
  //       AUTHOR         George Walton
  //       DATE WRITTEN   Extracted from AIRNET
  //       MODIFIED       Lixing Gu, 2/1/04
  //                      Revised the subroutine to meet E+ needs
  //       MODIFIED       Lixing Gu, 6/8/05
  //       RE-ENGINEERED  This subroutine is revised from CLVSKY developed by George Walton, NIST

  // PURPOSE OF THIS SUBROUTINE:
  // This subroutine solves simultaneous linear algebraic equations [A] * X = B
  // using L-U factored skyline form of [A] from "FACSKY"

  // METHODOLOGY EMPLOYED:
  // na

  // REFERENCES:
  // na

  // USE STATEMENTS:
  // na

  // Argument array dimensioning
  //AU.dim(IK(NetworkNumOfNodes + 1));
  //AD.dim(NetworkNumOfNodes);
  //AL.dim(IK(NetworkNumOfNodes + 1) - 1);
  //B.dim(NetworkNumOfNodes);
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

  int i;
  int JHK;
  int JHK1;
  int k;
  int LHK;
  Real64 SDOT;
  Real64 T1;

  // FLOW:
  JHK = 1;
  for (k = 2; k <= NEQ; ++k) {
    JHK1 = IK[k];
    LHK = JHK1 - JHK;
    if (LHK <= 0) continue;
    SDOT = 0.0;
    if (NSYM == 0) {
      for (i = 0; i <= LHK - 1; ++i) {
        SDOT += AU[JHK + i] * B[k - LHK + i];
      }
    } else {
      for (i = 0; i <= LHK - 1; ++i) {
        SDOT += AL[JHK + i] * B[k - LHK + i];
      }
    }
    B[k] -= SDOT;
    JHK = JHK1;
  }
  if (NSYM == 0) {
    for (k = 1; k <= NEQ; ++k) {
      B[k] *= AD[k];
    }
  }
  k = NEQ + 1;
  JHK1 = IK[k - 1];
  while (k != 1) {
    --k;
    if (NSYM == 1) B[k] *= AD[k];
    if (k == 1) break;
    //        IF(K.EQ.1) RETURN
    JHK = IK[k];
    T1 = B[k];
    for (i = 0; i <= JHK1 - JHK - 1; ++i) {
      B[k - JHK1 + JHK + i] -= AU[JHK + i] * T1;
    }
    JHK1 = JHK;
  }
}


void FACSKY(OneArray<Real64> AU,   // the upper triangle of [A] before and after factoring
  OneArray<Real64> AD,   // the main diagonal of [A] before and after factoring
  OneArray<Real64> AL,   // the lower triangle of [A] before and after factoring
  OneArray<int> const IK, // pointer to the top of column/row "K"
  int const NEQ,        // number of equations
  int const NSYM        // symmetry:  0 = symmetric matrix, 1 = non-symmetric
)
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
  int JHK;
  int JHK1;
  int LHK;
  int LHK1;
  int IMIN;
  int IMIN1;
  int JHJ;
  int JHJ1;
  int IC;
  int i;
  int j;
  int k;
  Real64 T1;
  Real64 T2;
  Real64 SDOT;
  Real64 SUMD;

  // FLOW:
  AD(1) = 1.0 / AD(1);
  JHK = 1;
  for (k = 2; k <= NEQ; ++k) {
    SUMD = 0.0;
    JHK1 = IK(k + 1);
    LHK = JHK1 - JHK;
    if (LHK > 0) {
      LHK1 = LHK - 1;
      IMIN = k - LHK1;
      IMIN1 = IMIN - 1;
      if (NSYM == 1) AL(JHK) *= AD(IMIN1);
      if (LHK1 != 0) {
        JHJ = IK(IMIN);
        if (NSYM == 0) {
          for (j = 1; j <= LHK1; ++j) {
            JHJ1 = IK(IMIN + j);
            IC = std::min(j, JHJ1 - JHJ);
            if (IC > 0) {
              SDOT = 0.0;
              for (i = 0; i <= IC - 1; ++i) {
                SDOT += AU(JHJ1 - IC + i) * AU(JHK + j - IC + i);
              }
              AU(JHK + j) -= SDOT;
            }
            JHJ = JHJ1;
          }
        } else {
          for (j = 1; j <= LHK1; ++j) {
            JHJ1 = IK(IMIN + j);
            IC = std::min(j, JHJ1 - JHJ);
            SDOT = 0.0;
            if (IC > 0) {
              for (i = 0; i <= IC - 1; ++i) {
                SDOT += AL(JHJ1 - IC + i) * AU(JHK + j - IC + i);
              }
              AU(JHK + j) -= SDOT;
              SDOT = 0.0;
              for (i = 0; i <= IC - 1; ++i) {
                SDOT += AU(JHJ1 - IC + i) * AL(JHK + j - IC + i);
              }
            }
            AL(JHK + j) = (AL(JHK + j) - SDOT) * AD(IMIN1 + j);
            JHJ = JHJ1;
          }
        }
      }
      if (NSYM == 0) {
        for (i = 0; i <= LHK1; ++i) {
          T1 = AU(JHK + i);
          T2 = T1 * AD(IMIN1 + i);
          AU(JHK + i) = T2;
          SUMD += T1 * T2;
        }
      } else {
        for (i = 0; i <= LHK1; ++i) {
          SUMD += AU(JHK + i) * AL(JHK + i);
        }
      }
    }
    if (AD(k) - SUMD == 0.0) {
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
    AD(k) = 1.0 / (AD(k) - SUMD);
    JHK = JHK1;
  }
}

void SLVSKY(OneArray<Real64> const AU, // the upper triangle of [A] before and after factoring
  OneArray<Real64> const AD, // the main diagonal of [A] before and after factoring
  OneArray<Real64> const AL, // the lower triangle of [A] before and after factoring
  OneArray<Real64> B,        // "B" vector (input); "X" vector (output).
  OneArray<int> const IK,     // pointer to the top of column/row "K"
  int const NEQ,            // number of equations
  int const NSYM            // symmetry:  0 = symmetric matrix, 1 = non-symmetric
)
{

  // SUBROUTINE INFORMATION:
  //       AUTHOR         George Walton
  //       DATE WRITTEN   Extracted from AIRNET
  //       MODIFIED       Lixing Gu, 2/1/04
  //                      Revised the subroutine to meet E+ needs
  //       MODIFIED       Lixing Gu, 6/8/05
  //       RE-ENGINEERED  This subroutine is revised from CLVSKY developed by George Walton, NIST

  // PURPOSE OF THIS SUBROUTINE:
  // This subroutine solves simultaneous linear algebraic equations [A] * X = B
  // using L-U factored skyline form of [A] from "FACSKY"

  // METHODOLOGY EMPLOYED:
  // na

  // REFERENCES:
  // na

  // USE STATEMENTS:
  // na

  // Argument array dimensioning
  //AU.dim(IK(NetworkNumOfNodes + 1));
  //AD.dim(NetworkNumOfNodes);
  //AL.dim(IK(NetworkNumOfNodes + 1) - 1);
  //B.dim(NetworkNumOfNodes);
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

  int i;
  int JHK;
  int JHK1;
  int k;
  int LHK;
  Real64 SDOT;
  Real64 T1;

  // FLOW:
  JHK = 1;
  for (k = 2; k <= NEQ; ++k) {
    JHK1 = IK(k + 1);
    LHK = JHK1 - JHK;
    if (LHK <= 0) continue;
    SDOT = 0.0;
    if (NSYM == 0) {
      for (i = 0; i <= LHK - 1; ++i) {
        SDOT += AU(JHK + i) * B(k - LHK + i);
      }
    } else {
      for (i = 0; i <= LHK - 1; ++i) {
        SDOT += AL(JHK + i) * B(k - LHK + i);
      }
    }
    B(k) -= SDOT;
    JHK = JHK1;
  }
  if (NSYM == 0) {
    for (k = 1; k <= NEQ; ++k) {
      B(k) *= AD(k);
    }
  }
  k = NEQ + 1;
  JHK1 = IK(k);
  while (k != 1) {
    --k;
    if (NSYM == 1) B(k) *= AD(k);
    if (k == 1) break;
    //        IF(K.EQ.1) RETURN
    JHK = IK(k);
    T1 = B(k);
    for (i = 0; i <= JHK1 - JHK - 1; ++i) {
      B(k - JHK1 + JHK + i) -= AU(JHK + i) * T1;
    }
    JHK1 = JHK;
  }
}

}
