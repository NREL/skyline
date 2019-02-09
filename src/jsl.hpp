// Copyright (c) 1996-2018, Jason W. DeGraw
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
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

} // jsl

#endif
