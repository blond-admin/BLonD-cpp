/*
 * openmp.h
 *
 *  Created on: Octomber 4, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDE_BLOND_OPENMP_H_
#define INCLUDE_BLOND_OPENMP_H_

#ifdef USE_OMP

#include <omp.h>

#else
namespace blond {
    static inline int omp_get_num_threads() {return 1;}
    static inline int omp_get_max_threads() {return 1;}
    static inline int omp_get_thread_num() {return 0;}
    static inline void omp_set_num_threads(int threads) {}
} // blond
#endif

#endif /* INCLUDE_BLOND_OPENMP_H_ */
