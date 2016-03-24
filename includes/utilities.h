/*
 * utilities.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_UTILITIES_H_
#define INCLUDES_UTILITIES_H_

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mm_malloc.h>
#include <sys/time.h>
#include "configuration.h"

#define TIMING
#define dprintf(...)    fprintf(stdout, __VA_ARGS__)     // Debug printf

#define CHECK_ERROR(a)                                       \
   if (a)                                                    \
   {                                                         \
      perror("Error at line\n\t" #a "\nSystem Msg");         \
      assert ((a) == 0);                                     \
   }

// sort an array with regards to another array
struct MyComparator {
	ftype *a;
	MyComparator(ftype* _a) :
			a(_a) {
	}

	bool operator()(int i1, int i2) {
		return a[i1] < a[i2];
	}
};

char const* GETENV(char const* envstr);

void *aligned_malloc(size_t n);

void delete_array(ftype *p);
void delete_array(bool *p);
void delete_array(int *p);

void dump(ftype* a, int n, const char* s);

void dump(int* a, int n, const char* s);

double time_diff(timespec const& end, timespec const& begin);

void get_time(timespec& ts);

timespec get_time();

double time_elapsed(timespec const& begin);

void print_time(char const* prompt, timespec const& begin, timespec const& end);

void print_time(char const* prompt, double diff);

void print_time_elapsed(char const* prompt, timespec const& begin);

#endif /* INCLUDES_UTILITIES_H_ */
