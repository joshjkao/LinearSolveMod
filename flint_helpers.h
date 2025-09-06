#pragma once
#define FLINT_INCLUDED

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>

// Wrappers for FLINT functions
// (mostly for debugging)

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
std::vector<std::vector<T>> FLINT_RREF_Modular(std::vector<std::vector<T>> &mat,
                                               const T &mod) {
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
