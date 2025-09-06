#pragma once
#include <algorithm>
#include <cmath>
#include <ranges>
#include <tuple>
#include <vector>

#ifdef DEBUG
#include "util.h"
#include <iostream>
#endif

#ifdef FLINT
#include "flint_helpers.h"
#endif

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

// Helper math functions
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
template <typename T> inline T fdiv(const T a, const T b) {
	return floor((double)a / b);
}
template <typename T> inline T PositiveMod(const T &a, const T &p) {
	T ret = a % p;
	if (ret < 0)
		ret += p;
	return ret;
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
template <typename T>
std::vector<T> MatMulMod(const std::vector<std::vector<T>> &mat,
                         const std::vector<T> &vec,
                         const std::vector<T> &moduli) {
	std::vector<T> ret(mat.size(), 0);
	for (size_t i = 0; i < mat.size(); ++i) {
		for (size_t j = 0; j < mat[0].size(); ++j) {
			ret[i] += mat[i][j] * vec[j];
		}
		if (moduli[i] != 0)
			ret[i] %= moduli[i];
	}
	return ret;
}
template <typename T>
std::vector<T> MatMulMod(const std::vector<std::vector<T>> &mat,
                         const std::vector<T> &vec, const T &mod) {
	std::vector<T> ret(mat.size(), 0);
	for (size_t i = 0; i < mat.size(); ++i) {
		for (size_t j = 0; j < mat[0].size(); ++j) {
			ret[i] += mat[i][j] * vec[j];
			ret[i] %= mod;
		}
	}
	return ret;
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
std::vector<std::vector<T>> RREF_Modular(std::vector<std::vector<T>> &A,
                                         const T &mod) {
	size_t m = A.size();
	if (m == 0)
		return A;
	size_t n = A[0].size();

	auto mat = A;

	size_t row = 0;
	for (size_t col = 0; col < n && row < m; ++col) {
		// identify pivot and swap if necessary
		// if no pivot, skip this column
		auto pivot =
		    std::ranges::find_if(mat | std::views::drop(row),
		                         [&](auto r) { return r[col] % mod != 0; });
		if (pivot == mat.end())
			continue;
		std::swap(*pivot, mat[row]);

		// reduce this row
		T inv = ModularInverse(mat[row][col] % mod, mod);
		if (inv == -1)
			continue;
		for (size_t j = col; j < n; ++j) {
			mat[row][j] = PositiveMod(mat[row][j] * inv, mod);
		}

		// eliminate above and below pivot
		for (size_t i = 0; i < m; ++i) {
			if (i == row)
				continue;
			T factor = mat[i][col];
			if (factor != 0) {
				for (size_t j = 0; j < n; ++j) {
					mat[i][j] =
					    PositiveMod(mat[i][j] - factor * mat[row][j], mod);
				}
			}
		}
		++row;
	}

	return mat;
}

// HNF Helper Functions
template <typename T>
void HNF_FixDiag(std::vector<T> &u, const T &a, const std::vector<T> &v,
                 const T &M, size_t m) {
	for (size_t i = 0; i < m; i++) {
		u[i] = (a * v[i]) % M;
	}
}
template <typename T>
void HNF_ReduceW(std::vector<T> &u, const T &a, const std::vector<T> &v,
                 const T &M, size_t m) {
	for (size_t i = 0; i < m; i++) {
		u[i] = (u[i] - a * v[i]) % M;
	}
}
template <typename T>
void HNF_EuclUpdate(std::vector<T> &u, std::vector<T> &v, const T &a,
                    const T &b, const T &c, const T &d, const T &M) {
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
				HNF_EuclUpdate(A[j], A[k], c1, c2, v, u, D);
			}
		}

		XGCD(d, u, v, A[k][i], D);
		HNF_FixDiag(W[i], u, A[k], D, i + 1);
		if (W[i][i] == 0)
			W[i][i] = D;

		for (j = i + 1; j < m; j++) {
			c1 = fdiv(W[j][i], W[i][i]);
			HNF_ReduceW(W[j], c1, W[i], D, i + 1);
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
			continue;
		}

		auto x = MatMulMod(H1, y, p);

		T M_new = M * p;
		for (auto &&[xi, Xi] : std::views::zip(x, X)) {
			T inv = ModularInverse(M, p);
			T delta = ((xi - Xi) * inv) % p;
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

template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli) {

	size_t num_zeros = std::ranges::count(moduli, 0);
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
		H = std::move(H1);
	}

	namespace rng = std::ranges;

	std::vector<T> soln;
	auto is_soln_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m),
		                   [](auto e) { return e == 0; }) &&
		       row[m] == 1;
	};
	auto soln_it = rng::find_if(H, is_soln_row);
	if (soln_it != H.end()) {
		for (const auto &e : *soln_it | rng::views::drop(m + 1)) {
			soln.push_back(e);
		}
	}

	std::vector<std::vector<T>> nulls;
	auto is_null_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m + 1),
		                   [](auto e) { return e == 0; });
	};
	auto null_rows = rng::views::filter(H, is_null_row);
	for (const auto &row : null_rows) {
		std::vector<T> null;
		for (const auto &e : row | rng::views::drop(m + 1)) {
			null.push_back(e);
		}
		nulls.push_back(null);
	}

	return {soln, nulls};
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
	size_t num_zeros = std::ranges::count(moduli, 0);
	auto nonzero_moduli =
	    moduli | std::views::filter([](auto e) { return e != 0; });

	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n;
	size_t augmat_n = m + n;
	size_t aug1_m = augmat_m;

	auto aug1 = std::vector<std::vector<T>>(aug1_m, std::vector<T>(aug1_m, 0));

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

	T d = 1;
	if (num_zeros > 0)
		d = Det(zero_block);
	for (const auto &m : nonzero_moduli)
		d *= m;

#ifdef DEBUG
	T d1 = d;
	for (const auto &m : moduli)
		d1 /= m;
	if (num_zeros > 0)
		d1 /= Det(zero_block);
#endif

	std::vector<std::vector<T>> H1, H;

	H1 = HNF_Modular(aug1, d);

	if (num_zeros > 0) {
		H = HNF_AddColumns(H1, augmat);
	} else {
		H = std::move(H1);
	}

	std::vector<std::vector<T>> nulls;

	namespace rng = std::ranges;
	auto is_null_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m),
		                   [](auto e) { return e == 0; });
	};
	auto null_rows = rng::views::filter(H, is_null_row);
	for (const auto &row : null_rows) {
		std::vector<T> null;
		for (const auto &e : row | rng::views::drop(m)) {
			null.push_back(e);
		}
		nulls.push_back(null);
	}

	return {nulls};
}

} // namespace LinSolveMod
