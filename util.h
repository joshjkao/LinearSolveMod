#pragma once
#include <iostream>
#include <vector>
#include <ranges>


template <typename T>
std::ostream& operator<<(std::ostream& os, 
    const std::vector<T>& vec)
{
    os << "[";
    for (auto& v: vec | std::views::take(vec.size()-1)) {
        os << v << ",";
    }
    os << vec.back() << "]";
    return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, 
    const std::vector<std::vector<T>>& mat)
{
    os << "[";
    for (auto& m: mat | std::views::take(mat.size()-1)) {
        os << m << "\n";
    }
    os << mat.back() << "]";
    return os;
}


template <typename T>
std::vector<T> operator*(const T& i, const std::vector<T>& vec) {
    auto ret = vec;
    for (auto& v: ret) v *= i;
    return ret;
}


template <typename T>
std::vector<T> operator+(
    const std::vector<T>& u, 
    const std::vector<T>& v) 
{
    auto ret = u;
    for (size_t i = 0; i < v.size(); i++) ret[i] += v[i];
    return ret;
}


template <typename T>
std::vector<T> MatMulMod(
    const std::vector<std::vector<T>>& mat, 
    const std::vector<T>& vec, 
    const std::vector<T>& moduli)
{
    std::vector<T> ret(mat.size(), 0);
    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[0].size(); ++j) {
            ret[i] += mat[i][j]*vec[j];
        }
        ret[i] %= moduli[i];
    }
    return ret;
}
