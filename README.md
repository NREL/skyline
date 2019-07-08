# skyline

[![Build Status](https://travis-ci.com/jasondegraw/skyline.svg?branch=master)](https://travis-ci.com/jasondegraw/skyline)

A header-only C++ direct solver for sparse matrices stored in skyline format. Include the header and then create the object:

```
std::vector<size_t> heights{ {0, 1, 1, 1, 1, 1, 1} };
skyline::SymmetricMatrix<size_t, float, std::vector> sky(heights);
```

Presently, only a LDLT decomposition solver is available:

```
sky.ldlt_solve(x);
```

The decomposition is done in place, so intermediate calls are available to solve multiple right hand sides. 
