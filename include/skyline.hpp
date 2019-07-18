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
#ifndef SKYLINE_HPP
#define SKYLINE_HPP

#include <numeric>
#include <optional>

namespace skyline {

template <typename I, typename R, template <typename ...> typename V> class SymmetricMatrix
{
public:

  SymmetricMatrix(V<V<R>> &M)
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

#ifdef SKYLINE_SINGLE_ARRAY
    // Copy into the am array
    m_am.resize(n + sum);
#else
    // Copy into the au/d arrays
    m_ad.resize(n);
    m_au.resize(sum);
#endif
    I count = 0;
    for (I i = 0; i < n; i++) {
#ifdef SKYLINE_SINGLE_ARRAY
      m_am[i] = M[i][i];
#else
      m_ad[i] = M[i][i];
#endif
      I j;
      for (j = 0; j < i; j++) {
        if (M[i][j] != 0.0) {
          break;
        }
      }
      for (; j < i; j++) {
#ifdef SKYLINE_SINGLE_ARRAY
        m_am[n + count] = M[i][j];
#else
        m_au[count] = M[i][j];
#endif
        ++count;
      }
    }

    m_v.resize(n);
    m_n = n;
  }

  SymmetricMatrix(V<I> &heights) : m_ih(heights)
  {
    I n = m_ih.size();
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

#ifdef SKYLINE_SINGLE_ARRAY
    // Size the am array
    I total{ n + sum };
    m_am.resize(total);
    for (I i = 0; i < total; i++) {
      m_am[i] = 0.0;
    }
#else
    // Size the au/d array
    m_ad.resize(n);
    for (I i = 0; i < n; i++) {
      m_ad[i] = 0.0;
    }

    m_au.resize(sum);
    for (I i = 0; i < sum; i++) {
      m_au[i] = 0.0;
    }
#endif

    m_v.resize(n);
    m_n = n;
  }

  void fill(R v = 0.0)
  {
#ifdef SKYLINE_SINGLE_ARRAY
    std::fill(m_am.begin(), m_am.end(), v);
#else
    std::fill(m_ad.begin(), m_ad.end(), v);
    std::fill(m_au.begin(), m_au.end(), v);
#endif
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
#ifdef SKYLINE_SINGLE_ARRAY
    return V<R>(m_am.begin(), m_am.begin() + m_n);
#else
    return m_ad;
#endif
  }

  V<R> upper() const
  {
#ifdef SKYLINE_SINGLE_ARRAY
    return V<R>(m_am.begin()+m_n, m_am.end());
#else
    return m_au;
#endif
  }

  V<R> lower() const
  {
#ifdef SKYLINE_SINGLE_ARRAY
    return V<R>(m_am.begin() + m_n, m_am.end());
#else
    return m_au;
#endif
  }

  R &operator()(I i)
  {
#ifdef SKYLINE_SINGLE_ARRAY
    return m_am[i];
#else
    return m_au[i];
#endif
  }

  R &diagonal(I i)
  {
#ifdef SKYLINE_SINGLE_ARRAY
    return m_am[i]; // This ends up being the same as operator()
#else
    return m_ad[i];
#endif
  }

  std::optional<I> index(I i, I j) const
  {
    if (m_im[j] <= i) {
      return m_n + m_ik[j] + i - m_im[j];
    } else if (m_im[j] == i) {
      return i;
    }
    return {};
  }

  virtual void utdu()
  {
    // j = 0, nothing much to do
    for (I k = 1; k < m_n; ++k) {
      if (m_im[k] == 0) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = m_n + m_ik[k];
        m_am[ij] /= m_am[0];
#else
        I ij = m_ik[k];
        m_au[ij] /= m_ad[0];
#endif
      }
    }
    // Now for the rest
    for (I j = 1; j < m_n; ++j) {
      // Compute v
      for (I i = 0; i < m_im[j]; ++i) {
        m_v[i] = 0.0;
      }
      for (I i = m_im[j]; i < j; ++i) {
#ifdef SKYLINE_SINGLE_ARRAY
        m_v[i] = m_am[m_n + m_ik[j] + i - m_im[j]] * m_am[i]; // OK, i >= m_im[j]
#else
        m_v[i] = m_au[m_ik[j] + i - m_im[j]] * m_ad[i]; // OK, i >= m_im[j]
#endif
      }
      // Compute the diagonal term
      R value = 0.0;
      for (I i = m_im[j]; i < j; ++i) {
#ifdef SKYLINE_SINGLE_ARRAY
        value += m_am[m_n + m_ik[j] + i - m_im[j]] * m_v[i];  // OK, i >= m_im[j]
#else
        value += m_au[m_ik[j] + i - m_im[j]] * m_v[i];  // OK, i >= m_im[j]
#endif
      }
#ifdef SKYLINE_SINGLE_ARRAY
      m_am[j] -= value;
#else
      m_ad[j] -= value;
#endif
      // Compute the rest of the row
      for (I k = j + 1; k < m_n; ++k) {
        if (m_im[k] <= j) {
          value = 0.0;
#ifdef SKYLINE_SINGLE_ARRAY
          for (I i = m_im[k]; i < j; ++i) {
            I ij = m_n + m_ik[k] + i - m_im[k]; // OK, i >= m_im[k]
            value += m_am[ij] * m_v[i];
          }
          I ij = m_n + m_ik[k] + j - m_im[k]; // OK, j >= m_im[k]
          m_am[ij] = (m_am[ij] - value) / m_am[j];
#else
          for (I i = m_im[k]; i < j; ++i) {
            I ij = m_ik[k] + i - m_im[k]; // OK, i >= m_im[k]
            value += m_au[ij] * m_v[i];
          }
          I ij = m_ik[k] + j - m_im[k]; // OK, j >= m_im[k]
          m_au[ij] = (m_au[ij] - value) / m_ad[j];
#endif
        }
      }
    }
  }

  virtual void forward_substitution(V<R> &b) const
  {
    // Solve Lz=b (Dy=z, Ux=y)
    for (I i = 1; i < m_n; ++i) {
      R value = 0.0;
      for (I k = m_im[i]; k < i; ++k) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = m_n + m_ik[i] + k - m_im[i];
        value += m_am[ij] * b[k];
#else
        I ij = m_ik[i] + k - m_im[i];
        value += m_au[ij] * b[k];
#endif
      }
      b[i] -= value;
    }
  }

  virtual void back_substitution(V<R> &z) const
  {
    // Account for the diagonal first (invert Dy=z)
    for (I j = 0; j < m_n; ++j) {
#ifdef SKYLINE_SINGLE_ARRAY
      z[j] /= m_am[j];
#else
      z[j] /= m_ad[j];
#endif
    }
    // Solve Ux=y
    for (I j = m_n - 1; j > 0; --j) {
      for (I k = m_im[j]; k < j; ++k) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = m_n + m_ik[j] + k - m_im[j];
        z[k] -= z[j] * m_am[ij];
#else
        I ij = m_ik[j] + k - m_im[j];
        z[k] -= z[j] * m_au[ij];
#endif
      }
    }
  }

  virtual void ldlt_solve(V<R> &b)
  {
    utdu();
    forward_substitution(b);
    back_substitution(b);
  }

  virtual void utdu_solve(V<R>& b)
  {
    utdu();
    forward_substitution(b);
    back_substitution(b);
  }

  I rows() const
  {
    return m_n;
  }

  I cols() const
  {
    return m_n;
  }

protected:

  I m_n;     // System size
  V<I> m_ik; // Index offsets to top of skylines
  V<I> m_ih; // Height of each skyline (not used, should probably be removed)
  V<I> m_im; // Minimum row, or top of skyline
#ifdef SKYLINE_SINGLE_ARRAY
  V<R> m_am; // The entire matrix in one vector, first the diagonal, then the rest
#else
  V<R> m_au; // Upper triangular part of matrix
  V<R> m_ad; // Diagonal of matrix
#endif
  V<R> m_v;  // Temporary used in solution
};

template <typename I, typename R, template <typename ...> typename V> class SymmetricSkipMatrix : public SymmetricMatrix<I, R, V>
{
public:

  SymmetricSkipMatrix(V<V<R>>& M) : SymmetricMatrix<I, R, V>(M)
  {
    m_skip.resize(this->m_n);
    m_ip.resize(this->m_n);
    for (I k = 0; k < this->m_n; ++k) {
      m_skip[k] = false;
      m_ip[k] = k;
    }
    m_locked = false;
    m_n_actual = this->m_n;
  }

  SymmetricSkipMatrix(V<I>& heights) : SymmetricMatrix<I, R, V>(heights)
  {
    m_skip.resize(this->m_n);
    m_ip.resize(this->m_n);
    I current = 0;
    for (I k = 0; k < this->m_n; ++k) {
      m_skip[k] = false;
      m_ip[k] = k;
    }
    m_locked = false;
    m_n_actual = this->m_n;
  }

  void utdu()
  {
    // j = 0, nothing much to do
    for (I k = 1; k < m_n_actual; ++k) {
      if (this->m_im[m_ip[k]] <= m_ip[0]) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = m_n + this->m_ik[m_ip[k]];
        this->m_am[ij] /= this->m_am[0];
#else
        I ij = this->m_ik[m_ip[k]];
        this->m_au[ij] /= this->m_ad[m_ip[0]];
#endif
      }
    }
    // Now for the rest
    for (I j = 1; j < m_n_actual; ++j) {
      // Compute v
      for (I i = 0; i < m_im[m_ip[j]]; ++i) {
        this->m_v[m_ip[i]] = 0.0;
      }
      for (I i = this->m_im[m_ip[j]]; i < j; ++i) {
#ifdef SKYLINE_SINGLE_ARRAY
        this->m_v[m_ip[i]] = this->m_am[m_n + this->m_ik[m_ip[j]] + i - this->m_im[m_ip[j]]] * this->m_am[m_ip[i]]; // OK, i >= m_im[j]
#else
        this->m_v[m_ip[i]] = this->m_au[this->m_ik[m_ip[j]] + i - this->m_im[m_ip[j]]] * this->m_ad[m_ip[i]]; // OK, i >= m_im[j]
#endif
      }
      // Compute the diagonal term
      R value = 0.0;
      for (I i = this->m_im[m_ip[j]]; i < j; ++i) {
#ifdef SKYLINE_SINGLE_ARRAY
        value += m_am[m_n + m_ik[m_ip[j]] + i - m_im[m_ip[j]]] * m_v[m_ip[i]];  // OK, i >= m_im[j]
#else
        value += this->m_au[this->m_ik[m_ip[j]] + i - this->m_im[m_ip[j]]] * m_v[this->m_ip[i]];  // OK, i >= m_im[j]
#endif
      }
#ifdef SKYLINE_SINGLE_ARRAY
      this->m_am[m_ip[j]] -= value;
#else
      this->m_ad[m_ip[j]] -= value;
#endif
      // Compute the rest of the row
      for (I k = j + 1; k < m_n_actual; ++k) {
        if (this->m_im[m_ip[k]] <= j) {
          value = 0.0;
#ifdef SKYLINE_SINGLE_ARRAY
          for (I i = this->m_im[k]; i < j; ++i) {
            I ij = this->m_n + this->m_ik[k] + i - this->m_im[k]; // OK, i >= m_im[k]
            value += this->m_am[ij] * this->m_v[i];
          }
          I ij = this->m_n + this->m_ik[k] + j - this->m_im[k]; // OK, j >= m_im[k]
          this->m_am[ij] = (this->m_am[ij] - value) / this->m_am[j];
#else
          for (I i = this->m_im[m_ip[k]]; i < j; ++i) {
            I ij = this->m_ik[m_ip[k]] + i - this->m_im[m_ip[k]]; // OK, i >= m_im[k]
            value += this->m_au[ij] * this->m_v[m_ip[i]];
          }
          I ij = this->m_ik[m_ip[k]] + j - this->m_im[m_ip[k]]; // OK, j >= m_im[k]
          this->m_au[ij] = (this->m_au[ij] - value) / this->m_ad[m_ip[j]];
#endif
        }
      }
    }
  }

  void forward_substitution(V<R>& b) const
  {
    // Solve Lz=b (Dy=z, Ux=y)
    for (I i = 1; i < m_n_actual; ++i) {
      R value = 0.0;
      for (I k = m_im[m_ip[i]]; k < i; ++k) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = this->m_n + this->m_ik[m_ip[i]] + k - this->m_im[m_ip[i]];
        value += this->m_am[ij] * b[m_ip[k]];
#else
        I ij = this->m_ik[m_ip[i]] + k - this->m_im[m_ip[i]];
        value += this->m_au[ij] * b[m_ip[k]];
#endif
      }
      b[m_ip[i]] -= value;
    }
  }

  void back_substitution(V<R>& z) const
  {
    // Account for the diagonal first (invert Dy=z)
    for (I j = 0; j < m_n_actual; ++j) {
#ifdef SKYLINE_SINGLE_ARRAY
      z[m_ip[j]] /= this->m_am[m_ip[j]];
#else
      z[m_ip[j]] /= this->m_ad[m_ip[j]];
#endif
    }
    // Solve Ux=y
    for (I j = m_n_actual - 1; j > 0; --j) {
      for (I k = this->m_im[m_ip[j]]; k < j; ++k) {
#ifdef SKYLINE_SINGLE_ARRAY
        I ij = this->m_n + this->m_ik[m_ip[j]] + k - this->m_im[m_ip[j]];
        z[m_ip[k]] -= z[m_ip[j]] * this->m_am[ij];
#else
        I ij = this->m_ik[m_ip[j]] + k - this->m_im[m_ip[j]];
        z[m_ip[k]] -= z[m_ip[j]] * this->m_au[ij];
#endif
      }
    }
  }

  virtual void ldlt_solve(V<R>& b)
  {
    lock();
    utdu();
    forward_substitution(b);
    back_substitution(b);
    unlock();
  }

  virtual void utdu_solve(V<R>& b)
  {
    lock();
    utdu();
    forward_substitution(b);
    back_substitution(b);
    unlock();
  }

  bool skip(I i)
  {
    if (m_locked) {
      return false;
    }
    m_skip[i] = !m_skip[i];
    return true;
  }

  bool &skip(I i) const
  {
    return m_skip[i];
  }

  V<bool> skip()
  {
    return m_skip;
  }

  void unskip()
  {
    for (I i = 0; i < m_n; ++i) {
      m_skip[i] = false;
    }
  }

  V<I> ip()
  {
    return m_ip;
  }

  void lock()
  {
    if (m_locked) {
      return;
    }
    I current = 0;
    for (I i = 0; i < m_n; ++i) {
      if (!m_skip[i]) {
        m_ip[current] = i;
        ++current;
      }
    }
    m_n_actual = current;
  }

  void unlock()
  {
    m_locked = false;
  }

private:
  bool m_locked;
  V<bool> m_skip;
  V<I> m_ip;
  I m_n_actual;
};


}

#endif // !SKYLINE_HPP
