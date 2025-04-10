# C++ Implementation of Linear Solve with Multiple Moduli

## Running the code

This project requires the [Fast Library for Number Theory (FLINT)](https://flintlib.org/). Optionally, the unit tests are run through [GoogleTest](https://github.com/google/googletest).

The example can be compiled and run with:
```
mkdir build
make run
```

The unit tests can be run with:
```
make gtests
./build/gtests
```

**Notes**: You may need to edit the Makefile depending on how you've built Googletest and FLINT. The Makefile included here works with Ubuntu 24.04. On OSX 14.3, I ran into issues with Clang and FLINT, and then again with gcc and the Homebrew version of Googletest. These were resolved by manually building Googletest. 

## Code Interface

The important code is contained in `linsolvemod.h`, which provides two functions:

```
template<typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> LinSolveMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& rhs,
    const std::vector<T>& moduli
);

template<typename T>
std::vector<std::vector<T>> NullSpaceMultiMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& moduli
);
```

The first one is the reimplementation of Mathematica's LinearSolveMod function. The second one is just a wrapper around the first.

## How it works

We desire a solution to 
```math
Ax=b
```
where the i-th row of A and b is modulo m_i. To solve this, we construct a matrix of the form:
```math
\begin{bmatrix}
(-b) && \\
A^T && I && \\
diag(m_i) && \\
\end{bmatrix}
```
and compute its row-style Hermite Normal Form. This makes it upper triangular, and the solutions can be read off from the first n+1 entries in each row (where n is the size of b). 

If the first n entries of a row are 0 and the (n+1)-th entry is 1, then the next n entries contain coefficients of the linear combination of the columns of A that add to b modulo m_i. 

Similarly, if all of the first (n+1) entries of a row are 0, then the next n entries contain coefficients of a linear combination of the columns of A that add to the zero vector modulo m_i.
