#include "linsolvemod.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <random>
#include <vector>

typedef long arrtype;
namespace LSM = LinSolveMod;

std::vector<arrtype> NullVector(size_t n) { return std::vector<arrtype>(n, 0); }

TEST(LINSOLVEMOD_BASIC, THREEBYTHREE) {
	std::vector<std::vector<arrtype>> mat = {{1, 1, 0}, {0, 1, 2}, {4, 1, 3}};
	std::vector<arrtype> moduli = {2, 2, 3};
	std::vector<arrtype> rhs = {0, 0, 1};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, THREEBYTHREE_LARGE) {
	std::vector<std::vector<arrtype>> mat = {{4, 2, 1}, {7, 1, 23}, {6, 2, 11}};
	std::vector<arrtype> moduli = {2, 2, 3};
	std::vector<arrtype> rhs = {0, 0, 1};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, THREEBYTHREE_LARGER) {
	std::vector<std::vector<arrtype>> mat = {
	    {23, 9, 123}, {54, 12, 97}, {45, 7, 12}};
	std::vector<arrtype> moduli = {45, 12, 94};
	std::vector<arrtype> rhs = {10, 8, 1};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, THREEBYTHREE_NEGATIVES) {
	std::vector<std::vector<arrtype>> mat = {
	    {5, -7, 23}, {34, 6, 2}, {-9, 12, -5}};
	std::vector<arrtype> moduli = {5, 7, 9};
	std::vector<arrtype> rhs = {4, 0, 3};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, FIVEBYFIVE) {
	std::vector<std::vector<arrtype>> mat = {{1, 1, 0, 5, 6},
	                                         {0, 1, 2, 9, 23},
	                                         {4, 1, 3, 3, 8},
	                                         {45, 2, 4, 5, 6},
	                                         {6, 2, 44, 7, 8}};
	std::vector<arrtype> rhs = {0, 0, 1, 0, 0};
	std::vector<arrtype> moduli = {2, 2, 3, 1, 4};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, SIXBYSIX) {
	std::vector<std::vector<arrtype>> mat = {
	    {1, 1, 0, 5, 6, 3},  {0, 1, 2, 9, 23, 5},  {4, 1, 3, 3, 8, 9},
	    {45, 2, 4, 5, 6, 1}, {6, 2, 44, 7, 8, 45}, {56, 2, 4, 6, 2, 3}};
	std::vector<arrtype> rhs = {0, 0, 1, 0, 0, 0};
	std::vector<arrtype> moduli = {2, 2, 3, 1, 4, 3};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(null.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_BASIC, EXTREME) {
	std::random_device rnd_device;
	std::mt19937 eng(rnd_device());
	std::uniform_int_distribution<arrtype> distr(-100, 100);
	std::uniform_int_distribution<arrtype> mod(2, 100);
	auto gen = [&]() { return distr(eng); };
	auto genmod = [&]() { return mod(eng); };
	std::vector<std::vector<arrtype>> mat(100, std::vector<arrtype>(100));
	for (auto &row : mat) {
		std::generate(row.begin(), row.end(), gen);
	}
	std::vector<arrtype> rhs(100);
	std::generate(rhs.begin(), rhs.end(), gen);
	std::vector<arrtype> moduli(100);
	std::generate(moduli.begin(), moduli.end(), genmod);
	for (auto &&[r, m] : std::views::zip(rhs, moduli)) {
		r %= m;
		if (r < 0)
			r += m;
	}
	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(mat.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_RECTANGULAR, FOURBYSIX) {
	std::vector<std::vector<arrtype>> mat = {{3, 3, 1, 4, 6, 7},
	                                         {0, 1, 0, 0, 4, 9},
	                                         {0, 0, 19, 16, 2, 43},
	                                         {0, 0, 0, 3, 7, 6}};
	std::vector<arrtype> rhs = {0, 1, 0, 1};
	std::vector<arrtype> moduli = {3, 5, 7, 10};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(mat.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_INFMOD, FOURBYSIX) {
	std::vector<std::vector<arrtype>> mat = {{3, 3, 1, 4, 6, 7},
	                                         {0, 1, 0, 0, 4, 9},
	                                         {0, 0, 19, 16, 2, 43},
	                                         {0, 0, 0, 3, 7, 6}};
	std::vector<arrtype> rhs = {0, 4, 6, 1};
	std::vector<arrtype> moduli = {3, 5, 7, 0};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(mat.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_INFMOD, SIXBYSIX) {
	std::vector<std::vector<arrtype>> mat = {
	    {7, 91, -17, 3, 11, 5}, {3, 34, 1, -5, 3, 7},    {-12, 4, 6, 13, 6, 7},
	    {-4, 7, 9, 11, -4, 6},  {-6, 4, 28, -12, 4, 67}, {7, -9, 23, 4, 52, 1}};
	std::vector<arrtype> rhs = {0, 1, 2, 3, 4, 5};
	std::vector<arrtype> moduli = {3, 5, 7, 0, 0, 0};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(mat.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

TEST(LINSOLVEMOD_INFMOD, SIXBYSIX_ZDET) {
	std::vector<std::vector<arrtype>> mat = {
	    {7, 91, -17, 3, 11, 5}, {3, 34, 1, -5, 3, 7},   {-12, 4, 6, 13, 6, 7},
	    {-4, 7, 9, 11, -4, 0},  {-6, 4, 28, -12, 4, 0}, {7, -9, 23, 4, 52, 0}};
	std::vector<arrtype> rhs = {0, 1, 2, 3, 4, 5};
	std::vector<arrtype> moduli = {3, 5, 7, 0, 0, 0};

	auto [soln, nulls] = LSM::LinSolveMod(mat, rhs, moduli);
	EXPECT_EQ(rhs, LSM::MatMulMod(mat, soln, moduli));
	for (const auto &null : nulls) {
		EXPECT_EQ(NullVector(mat.size()), LSM::MatMulMod(mat, null, moduli));
	}
}

int main(int argc, char **argv) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
