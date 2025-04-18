#include "linsolvemod.h"
#include "util.h"
#include <climits>


int main()
{
    std::vector<std::vector<long>> mat = {
        {1,1,0,5,6,3},
        {0,1,2,9,23,5},
        {4,1,3,3,8,9},
        {45,2,4,5,6,1},
        {6,2,44,7,8,45},
        {56,2,4,6,2,3}
    };
    std::vector<long> rhs = {0,0,0,0,0,0};
    std::vector<long> moduli = {2,2,3,1,4,3};

    mat = {{3,3,1,4,6,7},{0,1,0,0,4,9},{0,0,19,16,2,43},{0,0,0,3,7,6}};
    rhs = {0,1,0,1};
    moduli = {3,5,7,11};

    auto [soln, nulls] = LinSolveMod(mat, rhs, moduli);

    std::cout << "solution:\n" << soln << "\n";
    std::cout << "nulls:\n" << nulls << "\n\n";

    std::cout << "check the solution (should equal rhs):\n";
    std::cout << MatMulMod(mat, soln, moduli) << "\n";
    std::cout << "check the nulls (should equal zero):\n";
    for (const auto& null: nulls) {
        std::cout << MatMulMod(mat, null, moduli) << "\n";
    }
}