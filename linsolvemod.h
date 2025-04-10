#pragma once
#include <vector>
#include <tuple>
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"


template<typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> LinSolveMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& rhs,
    const std::vector<T>& moduli
)
{
    size_t m = mat.size();
    size_t n = mat[0].size();
    size_t augmat_m = n+moduli.size()+1;
    size_t augmat_n = m+n+moduli.size()+1;

    fmpz_mat_t A, H;
    fmpz_mat_init(A, augmat_m, augmat_n);
    fmpz_mat_init(H, augmat_m, augmat_n);

    // -rhs
    for (size_t j = 0; j < m; ++j) {
        fmpz_set_si(fmpz_mat_entry(A,0,j), -rhs[j]);
    }
    // row join with mat transpose
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            fmpz_set_si(fmpz_mat_entry(A,j+1,i), mat[i][j]);
        }
    }
    // row join with moduli diagonal
    for (size_t i = 0; i < moduli.size(); ++i) {
        fmpz_set_si(fmpz_mat_entry(A,1+n+i,i), moduli[i]);
    }
    // column join with identity
    for (size_t i = 0; i < augmat_m; ++i) {
        fmpz_set_si(fmpz_mat_entry(A,i,m+i), 1);
    }

#ifdef DEBUG
    fmpz_mat_print_pretty(A);
#endif

    fmpz_mat_hnf(H, A);

#ifdef DEBUG
    fmpz_mat_print_pretty(H);
#endif

    std::vector<T> soln;
    for (size_t i = 0; i < augmat_m; ++i) {
        bool isSoln = true;
        for (size_t j = 0; j < m; j++) {
            if (!fmpz_is_zero(fmpz_mat_entry(H,i,j))) isSoln = false;
        }
        if (isSoln && fmpz_is_one(fmpz_mat_entry(H,i,m))) {
            for (size_t j = m+1; j < m+n+1; ++j) {
                soln.push_back(fmpz_get_si(fmpz_mat_entry(H,i,j)));
            }
        }
    }

    std::vector<std::vector<T>> nulls;
    for (size_t i = 0; i < augmat_m; ++i) {
        bool isNull = true;
        for (size_t j = 0; j < m+1; ++j) {
            if (!fmpz_is_zero(fmpz_mat_entry(H,i,j))) isNull = false;
        }
        if (isNull) {
            std::vector<T> null;
            for (size_t j = m+1; j < m+n+1; ++j) {
                null.push_back(fmpz_get_si(fmpz_mat_entry(H,i,j)));
            }
            nulls.push_back(null);
        }
    }

    fmpz_mat_clear(A);
    fmpz_mat_clear(H);

    return {soln, nulls};

}


template<typename T>
std::vector<std::vector<T>> NullSpaceMultiMod(
    const std::vector<std::vector<T>>& mat,
    const std::vector<T>& moduli
)
{
    auto rhs = std::vector<T>(mat.size(), 0);
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    return nulls;
}