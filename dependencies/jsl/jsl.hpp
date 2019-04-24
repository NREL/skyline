// Copyright (c) 1996-2018, Jason W. DeGraw
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

#ifndef JSL_HPP
#define JSL_HPP

namespace jsl {

template <typename I, typename R, template <typename ...> typename V> bool is_symmetric(V<V<R>> &M)
{
  I N = M.size();
  for (I j = 0; j < N; ++j) {
    for (I i = j+1; i < N; ++i) {
      if (M[i][j] != M[j][i]) {
        //std::cout << i << ' ' << j << ' ' << M[i][j] << " != " << M[j][i] << std::endl;
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
template <typename I, typename R, template <typename ...> typename V>
  void GEnxn(I n, V<V<R>> &A, V<R> &x, V<R> &b, V<R> &z, V<I> &ip)
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

/*****************************************************************************/
/*                                                                           */
/*   LDLT - Compute the LDLT decomposition of a symmetric matrix.            */
/*                                                                           */
/*   Arguments:                                                              */
/*      int_t n ------------ number of rows and columns.                     */
/*      real_t A[][] ------- symmetric matrix.                               */
/*      real_t *v ---------- real work vector at least n in length.          */
/*                                                                           */
/*   This is a straightforward implementation of algorithm 4.1.2 from Golub  */
/*   and Van Loan:                                                           */
/*                                                                           */
/*   Golub and Van Loan, Matrix Computations (3rd Edition), The Johns        */
/*   Hopkins University Press, 1996, pg 139.                                 */
/*                                                                           */
/*****************************************************************************/
template <typename I, typename R, template <typename ...> typename V>
  void LDLT(I n, V<V<R>> &A, V<R> &v)
{
  // Algorithm:
  //
  // for j=0,n-1
  //   # Compute v[0:j]
  //   for i = 0:j-1
  //     v[i] = A[j][i]*A[i][i]
  //   end
  //   v[j] = A[j][j] - A[j][0:j-1]*v[0:j-1]
  //   # Store d[j] and compute L[j+1:n-1][j]
  //   A[j][j] = v[j]
  //   A[j+1:n-1][j] = (A[j+1:n-1][j] - A[j+1:n-1][0:j-1]*v[0:j-1])/v[j]
  // end
  //
  // j = 0, nothing much to do
  for (I k = 1; k < n; ++k) {
    A[k][0] /= A[0][0];
  }
  // Now for the rest
  for (I j = 1; j < n; ++j) {
    for (I i = 0; i < j; ++i) {
      v[i] = A[j][i] * A[i][i];
    }
    R value = 0.0;
    for (I i = 0; i < j; ++i) {
      value += A[j][i] * v[i];
    }
    v[j] = A[j][j] - value;
    A[j][j] = v[j];
    for (I k = j + 1; k < n; ++k) {
      value = 0.0;
      for (I i = 0; i < j; ++i) {
        value += A[k][i] * v[i];
      }
      A[k][j] = (A[k][j] - value) / v[j];
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*   UTDU - Compute the LDLT decomposition of a symmetric matrix in upper    */
/*          triangular form.                                                 */
/*                                                                           */
/*   Arguments:                                                              */
/*      int_t n ------------ number of rows and columns.                     */
/*      real_t A[][] ------- symmetric matrix.                               */
/*      real_t *v ---------- real work vector at least n in length.          */
/*                                                                           */
/*   This is a modified implementation of the algorithm from Golub           */
/*   and Van Loan that returns the upper triangular part:                    */
/*                                                                           */
/*   Golub and Van Loan, Matrix Computations (3rd Edition), The Johns        */
/*   Hopkins University Press, 1996, pg 139.                                 */
/*                                                                           */
/*****************************************************************************/
template <typename I, typename R, template <typename ...> typename V>
  void UTDU(I n, V<V<R>> &A, V<R> &v)
{
  // Algorithm:
  //
  // for j=0,n-1
  //   # Compute v[0:j]
  //   for i = 0:j-1
  //     v[i] = A[i][j]*A[i][i]
  //   end
  //   v[j] = A[j][j] - A[0:j-1][j]*v[0:j-1]
  //   # Store d[j] and compute U[j][j+1:n-1]
  //   A[j][j] = v[j]
  //   A[j][j+1:n-1] = (A[j][j+1:n-1] - A[0:j-1][j+1:n-1]*v[0:j-1])/v[j]
  // end
  //
  // j = 0, nothing much to do
  for (I k = 1; k < n; ++k) {
    A[0][k] /= A[0][0];
  }
  // Now for the rest
  for (I j = 1; j < n; ++j) {
    for (I i = 0; i < j; ++i) {
      v[i] = A[i][j] * A[i][i];
    }
    R value = 0.0;
    for (I i = 0; i < j; ++i) {
      value += A[i][j] * v[i];
    }
    v[j] = A[j][j] - value;
    A[j][j] = v[j];
    for (I k = j + 1; k < n; ++k) {
      value = 0.0;
      for (I i = 0; i < j; ++i) {
        value += A[i][k] * v[i];
      }
      A[j][k] = (A[j][k] - value) / v[j];
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*   forward_substitution - Perform foward substitution.                     */
/*                                                                           */
/*   Arguments:                                                              */
/*      int_t n ------------ number of rows and columns.                     */
/*      real_t L[][] ------- Lower triangular matrix with ones on the        */
/*                           diagonal.                                       */
/*      real_t b[] --------- real RHS vector at least n in length.           */
/*                                                                           */
/*   This is a straightforward implementation of algorithm 3.1.1 from Golub  */
/*   and Van Loan:                                                           */
/*                                                                           */
/*   Golub and Van Loan, Matrix Computations (3rd Edition), The Johns        */
/*   Hopkins University Press, 1996, pg 89.                                  */
/*                                                                           */
/*****************************************************************************/
template <typename I, typename R, template <typename ...> typename V>
  void forward_substitution(I n, V<V<R>> &L, V<R> &b)
{
  for (I i = 1; i < n; ++i) {
    for (I j = 0; j < i; ++j) {
      b[i] -= L[i][j] * b[j];
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*   back_substitution - Perform back substitution.                          */
/*                                                                           */
/*   Arguments:                                                              */
/*      int_t n ------------ number of rows and columns.                     */
/*      real_t U[][] ------- Upper triangular matrix.                        */
/*      real_t b[] --------- real RHS vector at least n in length.           */
/*                                                                           */
/*   This is a straightforward implementation of algorithm 3.1.2 from Golub  */
/*   and Van Loan:                                                           */
/*                                                                           */
/*   Golub and Van Loan, Matrix Computations (3rd Edition), The Johns        */
/*   Hopkins University Press, 1996, pg 89.                                  */
/*                                                                           */
/*****************************************************************************/
  template <typename I, typename R, template <typename ...> typename V>
  void back_substitution(I n, V<V<R>> &U, V<R> &b)
  {
    b[n - 1] /= U[n - 1][n - 1];
    I i = n - 1;
    do {
      --i;
      for (I j = i+1; j < n; ++j) {
        b[i] -= U[i][j] * b[j];
      }
      b[i] /= U[i][i];
    } while (i > 0);
  }

} // jsl

#endif
