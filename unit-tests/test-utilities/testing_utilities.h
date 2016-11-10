/*
 * testing_utilities.h
 *
 *  Created on: Oct 28, 2016
 *      Author: kiliakis
 */

#ifndef TESTING_UTILITIES_H_
#define TESTING_UTILITIES_H_

#include <blond/configuration.h>
#include <algorithm>
#include <string>
#include <cmath>

static inline void ASSERT_NEAR_LOOP(const f_vector_t &refV,
                                    const f_vector_t &realV,
                                    const std::string &varName,
                                    const double epsilon = 1e-8)
{
    auto max = *std::max_element(refV.begin(), refV.end());
    max = std::abs(max);
    ASSERT_EQ(refV.size(), realV.size());
    for (uint i = 0; i < refV.size(); ++i) {
        auto ref = refV[i];
        auto real = realV[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of " << varName << " failed on i " << i << "\n";
    }
}


static inline void ASSERT_DOUBLE_EQ_LOOP(const f_vector_t &refV,
        const f_vector_t &realV,
        const std::string &varName)
{
    ASSERT_EQ(refV.size(), realV.size());
    for (uint i = 0; i < refV.size(); ++i) {
        auto ref = refV[i];
        auto real = realV[i];
        ASSERT_DOUBLE_EQ(ref, real)
                << "Testing of " << varName << " failed on i " << i << "\n";
    }
}

static inline void ASSERT_EQ_LOOP(const int_vector_t &refV,
                                  const int_vector_t &realV,
                                  const std::string &varName)
{
    ASSERT_EQ(refV.size(), realV.size());
    for (uint i = 0; i < refV.size(); ++i) {
        auto ref = refV[i];
        auto real = realV[i];
        ASSERT_EQ(ref, real)
                << "Testing of " << varName << " failed on i " << i << "\n";
    }
}


#endif /* TESTING_UTILITIES_H_ */
