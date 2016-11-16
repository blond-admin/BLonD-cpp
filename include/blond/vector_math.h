/*
 * vector_math.h
 *
 *  Created on: Nov 2, 2016
 *      Author: kiliakis
 */

#pragma once

#include <algorithm>
#include <functional>
#include <cassert>
#include <vector>
#include <iostream>
#include <iterator>

//
// Useful vector operations
//



// apply a unary operator op at every element of vector a (op a)
template <typename T, typename F>
static inline std::vector<T> apply_f(const std::vector<T> &a, F op)
{
    // const size_t size = a.size();
    // std::vector<T> result(size);
    // for (size_t i = 0; i < size; i++)
    //     result[i] = op(a[i]);

    std::vector<T> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(), op);
    return result;
}

// apply a unary operator op at every element of vector a (op a)
template <typename T, typename F>
static inline void apply_f_in_place(std::vector<T> &a, F op)
{
    // const size_t size = a.size();
    // for (size_t i = 0; i < size; i++)
    //     a[i] = op(a[i]);

    std::transform(a.begin(), a.end(), a.begin(), op);
}

// apply a binary operator op at every pair of elements from a, b (a op b)
template <typename T, typename F>
static inline std::vector<T> apply_f(const std::vector<T> &a,
                                     const std::vector<T> &b, F op)
{
    // const size_t size = a.size();
    // std::vector<T> result(size);
    // for (size_t i = 0; i < size; i++)
    //     result[i] = op(a[i], b[i]);

    assert(a.size() == b.size());
    std::vector<T> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), op);
    return result;
}

// apply a binary operator op at every pair of elements from a, b (a op b)
template <typename T, typename F>
static inline void apply_f_in_place(std::vector<T> &a,
                                    const std::vector<T> &b, F op)
{
    // const size_t size = a.size();
    // for (size_t i = 0; i < size; i++)
    //     a[i] = op(a[i], b[i]);

    assert(a.size() == b.size());
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), op);
}


//
// Vector - Vector basic operations
//

template <typename T>
static inline std::vector<T> operator+(const std::vector<T> &a,
                                       const std::vector<T> &b)
{
    return apply_f(a, b, std::plus<T>());
}

template <typename T>
static inline std::vector<T> operator+=(std::vector<T> &a,
                                        const std::vector<T> &b)
{
    apply_f_in_place(a, b, std::plus<T>());
    return a;
}


template <typename T>
static inline std::vector<T> operator-(const std::vector<T> &a,
                                       const std::vector<T> &b)
{
    return apply_f(a, b, std::minus<T>());
}

template <typename T>
static inline std::vector<T> operator-=(std::vector<T> &a,
                                        const std::vector<T> &b)
{
    apply_f_in_place(a, b, std::minus<T>());
    return a;
}

template <typename T>
static inline std::vector<T> operator*(const std::vector<T> &a,
                                       const std::vector<T> &b)
{
    return apply_f(a, b, std::multiplies<T>());
}

template <typename T>
static inline std::vector<T> operator*=(std::vector<T> &a,
                                        const std::vector<T> &b)
{
    apply_f_in_place(a, b, std::multiplies<T>());
    return a;
}


template <typename T>
static inline std::vector<T> operator/(const std::vector<T> &a,
                                       const std::vector<T> &b)
{
    return apply_f(a, b, std::divides<T>());
}

template <typename T>
static inline std::vector<T> operator/=(std::vector<T> &a,
                                        const std::vector<T> &b)
{
    apply_f_in_place(a, b, std::divides<T>());
    return a;
}


//
// Vector - scalar basic operations
//


template <typename T, typename U>
static inline std::vector<T> operator+(const std::vector<T> &a, const U &b)
{
    return apply_f(a, std::bind2nd(std::plus<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator+(const U &b, const std::vector<T> &a)
{
    return apply_f(a, std::bind2nd(std::plus<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator+=(std::vector<T> &a, const U &b)
{
    apply_f_in_place(a, std::bind2nd(std::plus<T>(), b));
    return a;
}

template <typename T, typename U>
static inline std::vector<T> operator-(const std::vector<T> &a, const U &b)
{
    return apply_f(a, std::bind2nd(std::minus<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator-(const U &b, const std::vector<T> &a)
{
    return apply_f(a, std::bind1st(std::minus<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator-=(std::vector<T> &a, const U &b)
{
    apply_f_in_place(a, std::bind2nd(std::minus<T>(), b));
    return a;
}

template <typename T, typename U>
static inline std::vector<T> operator*(const std::vector<T> &a, const U &b)
{
    return apply_f(a, std::bind2nd(std::multiplies<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator*(const U &b, const std::vector<T> &a)
{
    return apply_f(a, std::bind2nd(std::multiplies<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator*=(std::vector<T> &a, const U &b)
{
    apply_f_in_place(a, std::bind2nd(std::multiplies<T>(), b));
    return a;
}

template <typename T, typename U>
static inline std::vector<T> operator/(const std::vector<T> &a, const U &b)
{
    return apply_f(a, std::bind2nd(std::divides<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator/(const U &b, const std::vector<T> &a)
{
    return apply_f(a, std::bind1st(std::divides<T>(), b));
}

template <typename T, typename U>
static inline std::vector<T> operator/=(std::vector<T> &a, const U &b)
{
    apply_f_in_place(a, std::bind2nd(std::divides<T>(), b));
    return a;
}

template <typename T>
static inline T pick(const std::vector<T> &a, const int i)
{
    assert(std::abs(i) < a.size());
    if (i >= 0)
        return a[i];
    else
        return a.rbegin()[-i - 1];
}

template <typename T>
static inline std::vector<T> pick(const std::vector<T> &a,
                                  const std::vector<int> &indexes)
{
    std::vector<T> res;
    res.reserve(indexes.size());
    for (const auto &i : indexes) {
        res.push_back(a[i]);
    }
    return res;
}

template <typename T>
static inline std::vector<T> pick(const std::vector<T> &a,
                                  const std::vector<bool> &indexes)
{
    std::vector<T> res;
    res.reserve(indexes.size());
    for (uint i = 0; i < indexes.size(); i++)
        if (indexes[i] == true)
            res.push_back(a[i]);
    return res;
}

template <typename T>
static inline std::vector<T> slice(const std::vector<T> &a,
                                   const int start = 0,
                                   int end = 0,
                                   const int step = 1)
{
    if (end == 0) end = a.size();
    else if (end < 0) end = a.size() + end;
    std::vector<T> result;
    result.reserve((end - start + step - 1) / step);
    for (int i = start; i < end; i += step)
        result.push_back(a[i]);
    return result;
}

template <typename T>
static inline std::vector<bool> operator==(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::equal_to<T>());
    return result;
}

template <typename T>
static inline std::vector<bool> operator!=(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::not_equal_to<T>());
    return result;
}

template <typename T>
static inline std::vector<bool> operator<(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::less<T>());
    return result;
}

template <typename T>
static inline std::vector<bool> operator<=(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::less_equal<T>());
    return result;
}


template <typename T>
static inline std::vector<bool> operator>(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::greater<T>());
    return result;
}

template <typename T>
static inline std::vector<bool> operator>=(const std::vector<T> &a,
        const std::vector<T> &b)
{
    assert(a.size() == b.size());
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   result.begin(), std::greater_equal<T>());
    return result;
}

template <typename T>
static inline std::vector<bool> operator==(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::equal_to<T>(), b));
    return result;
}

template <typename T>
static inline std::vector<bool> operator!=(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::not_equal_to<T>(), b));
    return result;
}


template <typename T>
static inline std::vector<bool> operator<(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::less<T>(), b));
    return result;
}


template <typename T>
static inline std::vector<bool> operator<=(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::less_equal<T>(), b));
    return result;
}


template <typename T>
static inline std::vector<bool> operator>(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::greater<T>(), b));
    return result;
}


template <typename T>
static inline std::vector<bool> operator>=(const std::vector<T> &a, const T &b)
{
    std::vector<bool> result(a.size());
    std::transform(a.begin(), a.end(), result.begin(),
                   std::bind2nd(std::greater_equal<T>(), b));
    return result;
}


template <typename T>
static inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &a)
{
    os << "[ ";
    std::copy(a.begin(), a.end(), std::ostream_iterator<T>(os, ", "));
    os << "\b\b ]\n";
    return os;
}

