#pragma once
#include <cmath>
#include <ranges>
#include <tuple>
#include <vector>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>

#include "util.h"
#include <iostream>
#include <print>

namespace LinSolveMod {

// Solves the integer system of linear equations mat*x = rhs
// modulo the values in "moduli". Returns a solution to the system
// and a list of vectors spanning the null space of mat.
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli);

// Solves the integer system modulo a single value using Gaussian
// elimination
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

template <typename T> T ModularInverse(T a, T m) {
	T d, s, t;
	XGCD(d, s, t, a, m);
	if (d != 1)
		return -1;
	else {
		return (s % m + m) % m;
	}
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
	flint_rand_init(rand);
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

	// todo: how to pick primes?
	std::vector<T> primes = {17, 19, 23, 41, 53, 97, 89, 83};

	std::vector<T> X(m, 0);
	T M = 1;

	for (const auto &p : primes) {
		auto [y, _] = LinSolveMod(A1, a, p);
		if (y.empty()) {
			std::cout << "[HNF_AddColumn] Prime " << p << " was singular\n";
			continue;
		}

		auto x = MatMulMod(H1, y, p);

		T M_new = M * p;
		for (auto &&[xi, Xi] : std::views::zip(x, X)) {
			T inv = ModularInverse(M, p);
			if (inv == -1)
				std::cout << "[HNF_AddColumn] ModularInverse doesn't exist\n";
			T delta = ((xi - Xi) * ModularInverse(M, p)) % p;
			Xi += M * delta;
		}
		M = M_new;
	}

	return X;
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

	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < m; ++j)
			H[i][j] = H1[i][j];

	return H;
}

template <typename T> T Det(const std::vector<std::vector<T>> &A) {
	int n = A.size();
	if (n == 1)
		return A[0][0];
	if (n == 2)
		return A[0][0] * A[1][1] - A[0][1] * A[1][0];
	T d = 0;
	for (int c = 0; c < n; c++) {
		auto m = std::vector<std::vector<T>>(n - 1, std::vector<T>(n - 1));
		for (int i = 1; i < n; i++)
			for (int j = 0, k = 0; j < n; j++)
				if (j != c)
					m[i - 1][k++] = A[i][j];
		d += (c % 2 ? -1 : 1) * A[0][c] * Det(m);
	}
	return d;
}

template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli) {

	size_t num_zeros = std::count(moduli.begin(), moduli.end(), 0);
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

	std::vector<std::vector<T>> zero_block;
	for (size_t i = 0; i < num_zeros; ++i) {
		std::vector<T> zero_block_row;
		for (size_t j = 0; j < num_zeros; ++j) {
			zero_block_row.push_back(
			    augmat[i + n + 1 - num_zeros][j + m - num_zeros]);
		}
		zero_block.push_back(zero_block_row);
	}

	if (num_zeros > 0) {
		std::cout << augmat << "\n" << zero_block << "\n";
	}

	T d = 1;
	if (num_zeros > 0)
		d = Det(zero_block);
	for (const auto &m : nonzero_moduli)
		d *= m;

#ifdef DEBUG
	T d1 = d;
	for (const auto &m : nonzero_moduli)
		d1 /= m;
	if (num_zeros > 0)
		d1 /= Det(zero_block);
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
std::vector<std::vector<T>> RREF_Modular(std::vector<std::vector<T>> &mat,
                                         const T &mod) {
	std::cout << mat << "\n";
	size_t m = mat.size();
	size_t n = mat[0].size();

	fmpz_mat_t a;
	fmpz_mat_init(a, m, n);

	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			fmpz_set_si(fmpz_mat_entry(a, i, j), mat[i][j]);
		}
	}

	fmpz_mod_ctx_t ctx;
	fmpz_mod_ctx_init_ui(ctx, mod);

	fmpz_mod_mat_t A;
	fmpz_mod_mat_init(A, m, n, ctx);

	fmpz_mod_mat_set_fmpz_mat(A, a, ctx);

	fmpz_mod_mat_t R;
	fmpz_mod_mat_init(R, m, n, ctx);

	fmpz_mod_mat_rref(R, A, ctx);

	fmpz_mod_mat_get_fmpz_mat(a, R, ctx);

	std::vector<std::vector<T>> ret(m, std::vector<T>(n, 0));

	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			ret[i][j] = fmpz_get_si(fmpz_mat_entry(a, i, j));
		}
	}

	fmpz_mat_clear(a);
	fmpz_mod_mat_clear(A, ctx);
	fmpz_mod_mat_clear(R, ctx);
	fmpz_mod_ctx_clear(ctx);

	return ret;
}

template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const T &modulus) {

	size_t m = mat.size();
	size_t n = mat[0].size();

	std::vector<std::vector<T>> augmat(m, std::vector<T>(n + 1, 0));
	for (size_t i = 0; i < m; ++i)
		for (size_t j = 0; j < n; ++j)
			augmat[i][j] = mat[i][j];
	for (size_t i = 0; i < m; ++i)
		augmat[i][n] = rhs[i];
	std::cout << augmat << "\n";

	auto rref = RREF_Modular(augmat, modulus);

	std::vector<T> ret(m, 0);

	for (size_t i = 0; i < m; ++i) {
		if (rref[i][i] != 1)
			return {{}, {}};
		ret[i] = rref[i][n];
	}

	return {ret, {}};
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

} // namespace LinSolveMod
