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
		//dprintf("a[0] = %.12lf\n", a[0]);
	}

	bool operator()(ftype i1, ftype i2) {
		//dprintf("i1,i2 = %d, %d\n", i1, i2);

		return i1 < i2; //a[i1] < a[i2];
	}
};

static inline char const* GETENV(char const* envstr) {
	char const* env = getenv(envstr);
	if (!env)
		return "0";
	else
		return env;
}

static inline void *aligned_malloc(size_t n) {
	return _mm_malloc(n, 64);
}

static inline void delete_array(ftype *p) {
	if (p != NULL)
		delete[] p;
}

static inline void delete_array(bool *p) {
	if (p != NULL)
		delete[] p;
}

static inline void delete_array(int *p) {
	if (p != NULL)
		delete[] p;
}

static inline void dump(ftype* a, int n, const char* s) {
	dprintf("%s", s);
	for (int i = 0; i < n; ++i) {
		dprintf("%.12lf\n", a[i]);
	}
	dprintf("\n");

}

static inline void dump(int* a, int n, const char* s) {
	dprintf("%s", s);
	for (int i = 0; i < n; ++i) {
		dprintf("%d\n", a[i]);
	}
	dprintf("\n");

}

static inline double time_diff(timespec const& end, timespec const& begin) {
#ifdef TIMING
	double result;

	result = end.tv_sec - begin.tv_sec;
	result += (end.tv_nsec - begin.tv_nsec) / (double) 1000000000;

	return result;
#else
	return 0;
#endif
}

static inline void get_time(timespec& ts) {
#ifdef TIMING
	// volatile long noskip;
#if _POSIX_TIMERS > 0
	clock_gettime(CLOCK_REALTIME, &ts);
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	ts.tv_sec = tv.tv_sec;
	ts.tv_nsec = tv.tv_usec * 1000;
#endif
	//noskip = ts.tv_nsec;
#endif
}

static inline timespec get_time() {
	timespec t;
#ifdef TIMING
	get_time(t);
#endif
	return t;
}

static inline double time_elapsed(timespec const& begin) {
#ifdef TIMING
	timespec now;
	get_time(now);
	return time_diff(now, begin);
#else
	return 0;
#endif
}

static inline void print_time(char const* prompt, timespec const& begin,
		timespec const& end) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, time_diff(end, begin));
#endif
}

static inline void print_time(char const* prompt, double diff) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, diff);
#endif
}

static inline void print_time_elapsed(char const* prompt,
		timespec const& begin) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, time_elapsed(begin));
#endif
}

#endif /* INCLUDES_UTILITIES_H_ */
