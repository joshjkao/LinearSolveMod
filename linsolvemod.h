#pragma once
#include <cmath>
#include <flint/flint.h>
#include <numeric>
#include <ranges>
#include <tuple>
#include <vector>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "util.h"
#include <iostream>
#include <print>

// Solves the integer system of linear equations mat*x = rhs
// modulo the values in "moduli". Returns a solution to the system
// and a list of vectors spanning the null space of mat.
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli);

// Solves the integer system modulo a single value
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const T &moduli);

template <typename T> inline T fdiv(const T a, const T b) {
	return floor((double)a / b);
}

// Euclid extended algorithm
template <typename T> void XGCD(T &d, T &s, T &t, T a, T b) {
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
	if (aneg)
		s = -s;
	if (bneg)
		t = -t;
	d = a1;
}

// u = a*v % M
template <typename T>
void FixDiag(std::vector<T> &u, const T &a, const std::vector<T> &v, const T &M,
             size_t m) {
	for (size_t i = 0; i < m; i++) {
		u[i] = (a * v[i]) % M;
	}
}

// u = (u - a*v) % M
template <typename T>
void ReduceW(std::vector<T> &u, const T &a, const std::vector<T> &v, const T &M,
             size_t m) {
	for (size_t i = 0; i < m; i++) {
		u[i] = (u[i] - a * v[i]) % M;
	}
}

template <typename T>
void EuclUpdate(std::vector<T> &u, std::vector<T> &v, const T &a, const T &b,
                const T &c, const T &d, const T &M) {
	size_t m = u.size();

	T M1 = M >> 1;

	T t1, t2, t3;

	for (size_t i = 0; i < m; i++) {
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
template <typename T>
std::vector<std::vector<T>> HNF_Modular(const std::vector<std::vector<T>> &A_in,
                                        const T &D_in) {
	std::vector<std::vector<T>> A = A_in;

	size_t n = A.size();
	size_t m = A[0].size();

	T D = D_in;

	std::vector<std::vector<T>> W(m, std::vector<T>(m, 0));

	size_t i, j, k;
	T d, u, v, c1, c2;

	for (i = 0; i < m; ++i) {
		for (j = 0; j < n / 2; ++j) {
			std::swap(A[i][j], A[i][m - j - 1]);
		}
	}

	k = n - 1;

	for (i = m - 1; i + 1 > 0; i--) {
		for (j = k - 1; j + 1 > 0; j--) {
			if (A[j][i] != 0) {
				XGCD(d, u, v, A[k][i], A[j][i]);
				c1 = fdiv(A[k][i], d);
				c2 = fdiv(A[j][i], d);
				c2 = -c2;
				EuclUpdate(A[j], A[k], c1, c2, v, u, D);
			}
		}

		XGCD(d, u, v, A[k][i], D);
		FixDiag(W[i], u, A[k], D, i + 1);
		if (W[i][i] == 0)
			W[i][i] = D;

		for (j = i + 1; j < m; j++) {
			c1 = fdiv(W[j][i], W[i][i]);
			ReduceW(W[j], c1, W[i], D, i + 1);
		}

		D = fdiv(D, d);
		k--;
	}

	// fix matrix orientation
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n / 2; ++j) {
			std::swap(W[i][j], W[i][m - j - 1]);
		}
	}
	for (i = 0; i < m / 2; ++i) {
		std::swap(W[i], W[m - i - 1]);
	}

	return W;
}

template <typename T>
std::vector<std::vector<T>>
FLINT_HNF_PernetStein(const std::vector<std::vector<T>> &A_in) {
	size_t m = A_in.size();
	size_t n = A_in[0].size();
	fmpz_mat_t A, H;
	fmpz_mat_init(A, m, n);
	fmpz_mat_init(H, m, n);
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j)
			fmpz_set_si(fmpz_mat_entry(A, i, j), A_in[i][j]);
	flint_rand_t rand;
	flint_randinit(rand);
	fmpz_mat_hnf_pernet_stein(H, A, rand);
	auto ret = std::vector<std::vector<T>>(m, std::vector<T>(n, 0));
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j)
			ret[i][j] = fmpz_get_si(fmpz_mat_entry(H, i, j));
	fmpz_mat_clear(A);
	fmpz_mat_clear(H);
	return ret;
}

template <typename T>
std::vector<T> HNF_AddColumn(const std::vector<std::vector<T>> &H1,
                             const std::vector<std::vector<T>> &A1,
                             const std::vector<T> &a) {
	size_t m = a.size();

	std::vector<T> primes = {11, 13, 57};
	auto big_mod =
	    std::accumulate(primes.begin(), primes.end(), 1, std::multiplies<T>());

	std::vector<T> crt_acc(m, 0);

	std::vector<std::vector<T>> solns(m, std::vector<T>(primes.size(), 0));

	for (const auto &[i, p] : primes | std::views::enumerate) {
		auto [y, _] = LinSolveMod(A1, a, p);
		if (y.empty()) {
			std::cout << "[HNF_AddColumn] Prime " << p << " was singular\n";
			continue;
		}
		auto x = MatMulMod(H1, y, p);
		for (const auto &[j, xj] : x | std::views::enumerate) {
			solns[j][i] = xj;
		}
	}

	std::vector<std::vector<T>> equivalences(primes.size(),
	                                         std::vector<T>(1, 1));

	for (const auto &[i, row] : solns | std::views::enumerate) {
		auto [xi, _] = LinSolveMod(equivalences, row, primes);
		crt_acc[i] = xi[0];
	}

	for (auto &xi : crt_acc) {
		if (xi > big_mod / 2) {
			xi -= big_mod;
		}
	}

	return crt_acc;
}

// Returns the HNF of A, given H1 is the HNF of the square nonsingular
// submatrix given by the first m columns of A.
template <typename T>
std::vector<std::vector<T>>
HNF_AddColumns(const std::vector<std::vector<T>> &H1,
               const std::vector<std::vector<T>> &A) {
	size_t m = H1.size();
	size_t n = A[0].size();

	auto A1 = std::vector<std::vector<T>>(m, std::vector<T>(m, 0));
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < m; ++j)
			A1[i][j] = A[i][j];
	auto cols = std::vector<std::vector<T>>(n - m, std::vector<T>(m, 0));
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n - m; ++j)
			cols[j][i] = A[i][j + m];

	auto H = std::vector<std::vector<T>>(m, std::vector<T>(n, 0));

	for (const auto &[j, col] : cols | std::views::enumerate) {
		auto x = HNF_AddColumn(H1, A1, col);
		for (size_t i = 0; i < m; ++i) {
			H[i][j + m] = x[i];
		}
	}
	std::cout << H << "\n";

	std::cout << H1 << "\n";

	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < m; ++j)
			H[i][j] = H1[i][j];
	std::cout << FLINT_HNF_PernetStein(A) << "\n";
	std::cout << H << "\n";
	return H;
}

template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli) {

	auto num_zeros = std::count(moduli.begin(), moduli.end(), 0);
	auto nonzero_moduli =
	    moduli | std::views::filter([](auto e) { return e != 0; });

	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n + 1 - num_zeros;
	size_t augmat_n = m + n + 1;
	size_t aug1_m = augmat_m;

	auto aug1 = std::vector<std::vector<T>>(aug1_m, std::vector<T>(aug1_m, 0));

	auto augmat =
	    std::vector<std::vector<T>>(augmat_m, std::vector<T>(augmat_n, 0));

	// -rhs
	for (size_t j = 0; j < m; ++j)
		augmat[0][j] = -rhs[j];
	// row join with mat transpose
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			augmat[j + 1][i] = mat[i][j];
		}
	}
	// row join with moduli diagonal
	for (size_t i = 0; i < moduli.size() - num_zeros; ++i)
		augmat[1 + n + i][i] = moduli[i];
	// column join with identity
	for (size_t i = 0; i < n + 1; ++i)
		augmat[i][m + i] = 1;
	// take a square nonsingular submatrix
	for (size_t i = 0; i < aug1_m; ++i)
		for (size_t j = 0; j < aug1_m; ++j)
			aug1[i][j] = augmat[i][j];

	T d = 1;
	for (const auto &m : nonzero_moduli)
		d *= m;

	T d1 = d;
	for (const auto &m : nonzero_moduli)
		d1 /= m;

#ifdef DEBUG
	if (d1 != 1) {
		std::cout << "[LinSolveMod] Possible Overflow!\n";
	}
#endif

	std::vector<std::vector<T>> H1, H;

	H1 = HNF_Modular(aug1, d);
	if (num_zeros > 0) {
		H = HNF_AddColumns(H1, augmat);
	} else {
		H = H1;
	}

	std::vector<T> soln;
	for (size_t i = 0; i < augmat_m; ++i) {
		bool isSoln = true;
		for (size_t j = 0; j < m; j++) {
			if (H[i][j] != 0)
				isSoln = false;
		}
		if (isSoln && H[i][m] == 1) {
			for (size_t j = m + 1; j < m + n + 1; ++j) {
				soln.push_back(H[i][j]);
			}
		}
	}

	std::vector<std::vector<T>> nulls;
	for (size_t i = 0; i < augmat_m; ++i) {
		bool isNull = true;
		for (size_t j = 0; j < m + 1; ++j) {
			if (H[i][j] != 0)
				isNull = false;
		}
		if (isNull) {
			std::vector<T> null;
			for (size_t j = m + 1; j < m + n + 1; ++j) {
				null.push_back(H[i][j]);
			}
			nulls.push_back(null);
		}
	}

	return {soln, nulls};
}

template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const T &modulus) {

	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n + 1;
	size_t augmat_n = m + n + 1;

	auto augmat =
	    std::vector<std::vector<T>>(augmat_m, std::vector<T>(augmat_n, 0));

	// -rhs
	for (size_t j = 0; j < m; ++j)
		augmat[0][j] = -rhs[j];
	// row join with mat transpose
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			augmat[j + 1][i] = mat[i][j];
		}
	}
	// row join with moduli diagonal
	for (size_t i = 0; i < m; ++i)
		augmat[1 + n + i][i] = modulus;
	// column join with identity
	for (size_t i = 0; i < n + 1; ++i)
		augmat[i][m + i] = 1;

	std::vector<std::vector<T>> H;

	T big_mod = 1;
	for (size_t i = 0; i < m; ++i)
		big_mod *= modulus;

	#ifdef DEBUG
	T mod1 = big_mod;
	for (size_t i = 0; i < m; ++i) {
		mod1 /= modulus;
	}
	if (mod1 != 1) {
		std::cout << "[LinSolveMod (single modulus)] Possible Overflow!\n";
	}
	#endif

	H = HNF_Modular(augmat, big_mod);
	
	std::vector<T> soln;
	for (size_t i = 0; i < augmat_m; ++i) {
		bool isSoln = true;
		for (size_t j = 0; j < m; j++) {
			if (H[i][j] != 0)
				isSoln = false;
		}
		if (isSoln && H[i][m] == 1) {
			for (size_t j = m + 1; j < m + n + 1; ++j) {
				soln.push_back(H[i][j]);
			}
		}
	}

	std::vector<std::vector<T>> nulls;
	for (size_t i = 0; i < augmat_m; ++i) {
		bool isNull = true;
		for (size_t j = 0; j < m + 1; ++j) {
			if (H[i][j] != 0)
				isNull = false;
		}
		if (isNull) {
			std::vector<T> null;
			for (size_t j = m + 1; j < m + n + 1; ++j) {
				null.push_back(H[i][j]);
			}
			nulls.push_back(null);
		}
	}

	return {soln, nulls};
}

// Returns a list of vectors spanning the null space of mat.
template <typename T>
std::vector<std::vector<T>>
NullSpaceMultiMod(const std::vector<std::vector<T>> &mat,
                  const std::vector<T> &moduli) {
	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n;
	size_t augmat_n = m + n;

	auto augmat =
	    std::vector<std::vector<T>>(augmat_m, std::vector<T>(augmat_n, 0));

	// row join with mat transpose
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			augmat[j][i] = mat[i][j];
		}
	}
	// row join with moduli diagonal
	for (size_t i = 0; i < moduli.size(); ++i) {
		augmat[n + i][i] = moduli[i];
	}
	// column join with identity
	for (size_t i = 0; i < n; ++i) {
		augmat[i][m + i] = 1;
	}

	T d = 1;
	for (const auto &m : moduli)
		d *= m;

	T d1 = d;
	for (const auto &m : moduli)
		d1 /= m;

	std::vector<std::vector<T>> H;

	if (d == 0) {
		// rectangular matrix: use pernet stein
		H = FLINT_HNF_PernetStein(augmat);
	} else if (d1 != 1) {
		// possible overflow, use flint types
		H = FLINT_HNF_Modular(augmat, d);
	} else {
		// use native types
		H = HNF_Modular(augmat, d);
	}

	std::vector<std::vector<T>> nulls;
	for (size_t i = 0; i < augmat_m; ++i) {
		bool isNull = true;
		for (size_t j = 0; j < m; ++j) {
			if (H[i][j] != 0)
				isNull = false;
		}
		if (isNull) {
			std::vector<T> null;
			for (size_t j = m; j < m + n; ++j) {
				null.push_back(H[i][j]);
			}
			nulls.push_back(null);
		}
	}

	return {nulls};
}
