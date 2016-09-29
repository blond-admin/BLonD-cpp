/*
 * constants.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDE_BLOND_OMP_H_
#define INCLUDE_BLOND_OMP_H_

#ifdef USE_OMP

#include <omp.h>

#else

int omp_get_num_threads() {return 1;}
int omp_get_thread_num() {return 0;}
void omp_set_num_threads(int threads) {}

#endif

#endif /* INCLUDE_BLOND_OMP_H_ */
