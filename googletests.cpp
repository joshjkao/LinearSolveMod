#include <iostream>
#include <vector>
#include <ranges>
#include "linsolvemod.h"
#include "util.h"
#include "gtest/gtest.h"

typedef long arrtype;

std::vector<arrtype> NullVector(size_t n) {
    return std::vector<arrtype>(n, 0);
}

TEST(LINSOLVEMOD, THREEBYTHREE) {
    std::vector<std::vector<arrtype>> mat = {{1,1,0},{0,1,2},{4,1,3}};
    std::vector<arrtype> moduli = {2,2,3};
    std::vector<arrtype> rhs = {0,0,1};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

TEST(LINSOLVEMOD, THREEBYTHREE_LARGE) {
    std::vector<std::vector<arrtype>> mat = {{4,2,1},{7,1,23},{6,2,11}};
    std::vector<arrtype> moduli = {2,2,3};
    std::vector<arrtype> rhs = {0,0,1};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

TEST(LINSOLVEMOD, THREEBYTHREE_LARGER) {
    std::vector<std::vector<arrtype>> mat = {{23,9,123},{54,12,97},{45,7,12}};
    std::vector<arrtype> moduli = {45,12,94};
    std::vector<arrtype> rhs = {10,8,1};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

TEST(LINSOLVEMOD, THREEBYTHREE_NEGATIVES) {
    std::vector<std::vector<arrtype>> mat = {{5,-7,23},{34,6,2},{-9,12,-5}};
    std::vector<arrtype> moduli = {5,7,9};
    std::vector<arrtype> rhs = {4,0,3};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

TEST(LINSOLVEMOD, FIVEBYFIVE) {
    std::vector<std::vector<arrtype>> mat = {{1,1,0,5,6},{0,1,2,9,23},{4,1,3,3,8},{45,2,4,5,6},{6,2,44,7,8}};
    std::vector<arrtype> rhs = {0,0,1,0,0};
    std::vector<arrtype> moduli = {2,2,3,1,4};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

TEST(LINSOLVEMOD, SIXBYSIX) {
    std::vector<std::vector<arrtype>> mat = {{1,1,0,5,6,3},{0,1,2,9,23,5},{4,1,3,3,8,9},{45,2,4,5,6,1},{6,2,44,7,8,45},{56,2,4,6,2,3}};
    std::vector<arrtype> rhs = {0,0,1,0,0,0};
    std::vector<arrtype> moduli = {2,2,3,1,4,3};
    
    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);
    EXPECT_EQ(rhs, MatMulMod(mat, soln, moduli));
    for (const auto& null: nulls) {
        EXPECT_EQ(NullVector(null.size()), MatMulMod(mat, null, moduli));
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}