#pragma once
#include <algorithm>
#include <cmath>
#include <ranges>
#include <vector>

#ifdef DEBUG
#include "util.h"
#include <iostream>
#endif

#include "flint_helpers.h"

namespace LinSolveMod {

// Solves the integer system of linear equations mat*x = rhs
// modulo the values in "moduli". Returns a solution to the system
// and a list of vectors spanning the null space of mat.
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli);

// Returns a list of vectors spanning the null space of mat,
// whos columns are vectors are defined modulo the values in "moduli".
template <typename T>
std::vector<std::vector<T>>
NullSpaceMultiMod(const std::vector<std::vector<T>> &mat,
                  const std::vector<T> &moduli);


template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
LinSolveMod(const std::vector<std::vector<T>> &mat, const std::vector<T> &rhs,
            const std::vector<T> &moduli) {

	size_t num_zeros = std::ranges::count(moduli, 0);
	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n + 1 - num_zeros;
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
	for (size_t i = 0; i < moduli.size() - num_zeros; ++i)
		augmat[1 + n + i][i] = moduli[i];
	// column join with identity
	for (size_t i = 0; i < n + 1; ++i)
		augmat[i][m + i] = 1;

	auto H = FLINT_HNF_PernetStein(augmat);

	namespace rng = std::ranges;

	auto is_soln_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m),
		                   [](auto e) { return e == 0; }) &&
		       row[m] == 1;
	};
	auto soln_it = rng::find_if(H, is_soln_row);
	std::vector<T> soln;
	if (soln_it != H.end()) {
		soln = *soln_it | rng::views::drop(m + 1) | rng::to<std::vector<T>>();
	}

	auto is_null_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m + 1),
		                   [](auto e) { return e == 0; });
	};
	auto null_rows = rng::views::filter(H, is_null_row);
	std::vector<std::vector<T>> nulls;
	for (const auto &row : null_rows) {
		nulls.emplace_back(row | rng::views::drop(m + 1) |
		                   rng::to<std::vector<T>>());
	}

	return {soln, nulls};
}

template <typename T>
std::vector<std::vector<T>>
NullSpaceMultiMod(const std::vector<std::vector<T>> &mat,
                  const std::vector<T> &moduli) {
	size_t num_zeros = std::ranges::count(moduli, 0);
	size_t m = mat.size();
	size_t n = mat[0].size();
	size_t augmat_m = m + n - num_zeros;
	size_t augmat_n = m + n;

	std::vector<std::vector<T>> augmat(augmat_m, std::vector<T>(augmat_n, 0));

	// row join with mat transpose
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			augmat[j][i] = mat[i][j];
		}
	}
	// row join with moduli diagonal
	for (size_t i = 0; i < moduli.size() - num_zeros; ++i) {
		augmat[n + i][i] = moduli[i];
	}
	// column join with identity
	for (size_t i = 0; i < n; ++i) {
		augmat[i][m + i] = 1;
	}

	auto H = FLINT_HNF_PernetStein(augmat);

	namespace rng = std::ranges;
	auto is_null_row = [&](const auto &row) {
		return rng::all_of(row | rng::views::take(m),
		                   [](auto e) { return e == 0; });
	};
	auto null_rows = rng::views::filter(H, is_null_row);
	std::vector<std::vector<T>> nulls;
	for (const auto &row : null_rows) {
		nulls.emplace_back(row | rng::views::drop(m) |
		                   rng::to<std::vector<T>>());
	}

	return {nulls};
}

} // namespace LinSolveMod
