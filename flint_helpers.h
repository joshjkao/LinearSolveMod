#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <vector>

// Wrappers for FLINT functions

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
	flint_rand_clear(rand);
	return ret;
}
