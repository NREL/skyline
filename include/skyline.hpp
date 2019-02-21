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

    // Copy into the au/d array
    m_ad.resize(n);
    for (I i = 0; i < n; i++) {
      m_ad[i] = 0.0;
    }

    m_au.resize(sum);
    for (I i = 0; i < sum; i++) {
      m_au[i] = 0.0;
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

  R &operator()(I i)
  {
    return m_au[i];
  }

  R &diagonal(I i)
  {
    return m_ad[i];
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
        m_v[i] = m_au[m_ik[j] + i - m_im[j]] * m_ad[i]; // OK, i >= m_im[j]
      }
      // Compute the diagonal term
      R value = 0.0;
      for (I i = m_im[j]; i < j; ++i) {
        value += m_au[m_ik[j] + i - m_im[j]] * m_v[i];  // OK, i >= m_im[j]
      }
      m_ad[j] -= value;
      // Compute the rest of the row
      for (I k = j + 1; k < m_n; ++k) {
        if (m_im[k] <= j) {
          value = 0.0;
          for (I i = m_im[k]; i < j; ++i) {
            I ij = m_ik[k] + i - m_im[k]; // OK, i >= m_im[k]
            value += m_au[ij] * m_v[i];
          }
          I ij = m_ik[k] + j - m_im[k]; // OK, j >= m_im[k]
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

  I rows() const
  {
    return m_n;
  }

  I cols() const
  {
    return m_n;
  }

private:

  I m_n;     // System size
  V<I> m_ik; // Index offsets to top of skylines
  V<I> m_ih; // Height of each skyline (not used, should probably be removed
  V<I> m_im; // Minimum row, or top of skyline
  V<R> m_au; // Upper triangular part of matrix
  V<R> m_ad; // Diagonal of matrix
  V<R> m_v;  // Temporary used in solution
};

}

#endif // !SKYLINE_HPP
