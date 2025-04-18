#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#ifdef DEBUG
#include <iostream>
#include "util.h"
#endif


// Solves the integer system of linear equations mat*x = rhs
// modulo the values in "moduli". Returns a solution to the system
// and a list of vectors spanning the null space of mat.
template<typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> LinSolveMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& rhs,
    const std::vector<T>& moduli
);


// Returns a list of vector spanning the null space of mat.
template<typename T>
std::vector<std::vector<T>> NullSpaceMultiMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& moduli
);


template<typename T>
inline T fdiv(const T a, const T b) { return floor((double)a / b); }


// Euclid extended algorithm
template<typename T>
void XGCD(T& d, T& s, T& t, T a, T b) {
    s = 1, t = 0;
    bool aneg = false, bneg = false;
    if (a < 0) {
        a = -a;
        aneg = true;
    }
    if (b < 0) {
        b = -b;
        bneg = true;
    }
    T s1 = 0, t1 = 1, a1 = a, b1 = b;
    while (b1) {
        T q = a1 / b1;
        std::tie(s, s1) = std::make_tuple(s1, s - q * s1);
        std::tie(t, t1) = std::make_tuple(t1, t - q * t1);
        std::tie(a1, b1) = std::make_tuple(b1, a1 - q * b1);
    }
    if (aneg) s = -s;
    if (bneg) t = -t;
    d = a1;
}


template<typename T>
void FixDiag(
    std::vector<T>& u, 
    const T& a, 
    const std::vector<T>& v, 
    const T& M, 
    long m
)
{
    long i;
    T t1;

    for (i = 0; i < m; i++) {
        t1 = a * v[i];
        u[i] = t1 % M;
    }
}


template<typename T>
void ReduceW(
    std::vector<T>& u, 
    const T& a, 
    const std::vector<T>& v, 
    const T& M, long m
)
{
    long i;
    T t1, t2;

    for (i = 0; i < m; i++) {
        t1 = a * v[i];
        t2 = u[i] - t1;
        u[i] = t2 % M;
    }
}


template<typename T>
void EuclUpdate(
    std::vector<T>& u, 
    std::vector<T>& v, 
    const T& a, const T& b, 
    const T& c, const T& d, 
    const T& M
)
{
    long m = u.size(); 
    long i;

    T M1;
    M1 = M >> 1;

    T t1, t2, t3;

    for (i = 0; i < m; i++) {
        t1 = u[i] * a;
        t2 = v[i] * b;
        t1 += t2;
        t1 %= M;
        if (t1 > M1)
            t1 -= M;

        t3 = t1;

        t1 = u[i] * c;
        t2 = v[i] * d;
        t1 += t2;
        t1 %= M;
        if (t1 > M1)
            t1 -= M;

        u[i] = t3;
        v[i] = t1;
    }
}


// Compute the row-style Hermite Normal Form of A_in, where 
// D_in is the determinant of the lattice spanned by A_in
// This code adapted from NTL's implementation
template<typename T>
std::vector<std::vector<T>> HNF_Modular(
    const std::vector<std::vector<T>>& A_in, 
    const T& D_in
)
{
    std::vector<std::vector<T>> A = A_in;

    long n = A.size();
    long m = A[0].size();

    T D = D_in;

    auto W = std::vector<std::vector<T>>(m, std::vector<T>(m, 0));

    long i, j, k;
    T d, u, v, c1, c2;

    for (i = 0; i < m; ++i) {
        for (j = 0; j < n/2; ++j) {
            std::swap(A[i][j], A[i][m-j-1]);
        }
    }

    k = n-1;

    for (i = m-1; i >= 0; i--) {
        for (j = k-1; j >= 0; j--) {
            if (A[j][i] != 0) {
                XGCD(d, u, v, A[k][i], A[j][i]);
                c1 = fdiv(A[k][i], d);
                c2 = fdiv(A[j][i], d);
                c2 = -c2;
                EuclUpdate(A[j], A[k], c1, c2, v, u, D);
            }
        }

        XGCD(d, u, v, A[k][i], D);
        FixDiag(W[i], u, A[k], D, i+1);
        if (W[i][i] == 0) W[i][i] = D;

        for (j = i+1; j < m; j++) {
            c1 = fdiv(W[j][i], W[i][i]);
            ReduceW(W[j], c1, W[i], D, i+1);
        }

        D = fdiv(D, d);
        k--;
    }

    // fix matrix orientation
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n/2; ++j) {
            std::swap(W[i][j], W[i][m-j-1]);
        }
    }
    for (i = 0; i < m/2; ++i) {
        std::swap(W[i], W[m-i-1]);
    }

    return W;
}


template<typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> LinSolveMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& rhs,
    const std::vector<T>& moduli
)
{
    size_t m = mat.size();
    size_t n = mat[0].size();
    size_t augmat_m = m+n+1;
    size_t augmat_n = m+n+1;

    auto augmat = std::vector<std::vector<T>>(augmat_m, std::vector<T>(augmat_n, 0));

    // -rhs
    for (size_t j = 0; j < m; ++j) {
        augmat[0][j] = -rhs[j];
    }
    // row join with mat transpose
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmat[j+1][i] = mat[i][j];
        }
    }
    // row join with moduli diagonal
    for (size_t i = 0; i < moduli.size(); ++i) {
        augmat[1+n+i][i] = moduli[i];
    }
    // column join with identity
    for (size_t i = 0; i < n+1; ++i) {
        augmat[i][m+i] = 1;
    }

    T d = 1;
    for (const auto& m: moduli) d *= m;

#ifdef DEBUG
    T d1 = d;
    for (const auto& m: moduli) d1 /= m;
    if (d1 != 1) {
        std::cout << "WARNING: POSSIBLE OVERFLOW\n";
    }
    std::cout << "LinSolveMod: original augmat:\n";
    std::cout << augmat << "\n";
#endif

    auto H = HNF_Modular(augmat, d);

#ifdef DEBUG
    std::cout << "hnf:\n";
    std::cout << H << "\n";
#endif

    std::vector<T> soln;
    for (size_t i = 0; i < augmat_m; ++i) {
        bool isSoln = true;
        for (size_t j = 0; j < m; j++) {
            if (H[i][j] != 0) isSoln = false;
        }
        if (isSoln && H[i][m] == 1) {
            for (size_t j = m+1; j < m+n+1; ++j) {
                soln.push_back(H[i][j]);
            }
        }
    }

    std::vector<std::vector<T>> nulls;
    for (size_t i = 0; i < augmat_m; ++i) {
        bool isNull = true;
        for (size_t j = 0; j < m+1; ++j) {
            if (H[i][j] != 0) isNull = false;
        }
        if (isNull) {
            std::vector<T> null;
            for (size_t j = m+1; j < m+n+1; ++j) {
                null.push_back(H[i][j]);
            }
            nulls.push_back(null);
        }
    }

    return {soln, nulls};

}


template<typename T>
std::vector<std::vector<T>> NullSpaceMultiMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& moduli
)
{
    size_t m = mat.size();
    size_t n = mat[0].size();
    size_t augmat_m = m+n;
    size_t augmat_n = m+n;

    auto augmat = std::vector<std::vector<T>>(augmat_m, std::vector<T>(augmat_n, 0));

    // row join with mat transpose
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmat[j][i] = mat[i][j];
        }
    }
    // row join with moduli diagonal
    for (size_t i = 0; i < moduli.size(); ++i) {
        augmat[n+i][i] = moduli[i];
    }
    // column join with identity
    for (size_t i = 0; i < n; ++i) {
        augmat[i][m+i] = 1;
    }

    T d = 1;
    for (const auto& m: moduli) d *= m;

#ifdef DEBUG
    std::cout << "LinSolveMod: original augmat:\n";
    std::cout << augmat << "\n";
#endif

    auto H = HNF_Modular(augmat, d);

#ifdef DEBUG
    std::cout << "hnf:\n";
    std::cout << H << "\n";
#endif

    std::vector<std::vector<T>> nulls;
    for (size_t i = 0; i < augmat_m; ++i) {
        bool isNull = true;
        for (size_t j = 0; j < m; ++j) {
            if (H[i][j] != 0) isNull = false;
        }
        if (isNull) {
            std::vector<T> null;
            for (size_t j = m; j < m+n; ++j) {
                null.push_back(H[i][j]);
            }
            nulls.push_back(null);
        }
    }

    return {nulls};
}