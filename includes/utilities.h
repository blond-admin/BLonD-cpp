/*
 * utilities.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mm_malloc.h>
#include <sys/time.h>
#include "configuration.h"

#ifndef INCLUDES_UTILITIES_H_
#define INCLUDES_UTILITIES_H_
#define TIMING
#define dprintf(...)    fprintf(stdout, __VA_ARGS__)     // Debug printf

#define CHECK_ERROR(a)                                       \
   if (a)                                                    \
   {                                                         \
      perror("Error at line\n\t" #a "\nSystem Msg");         \
      assert ((a) == 0);                                     \
   }

static inline char const* GETENV(char const* envstr) {
	char const* env = getenv(envstr);
	if (!env)
		return "0";
	else
		return env;
}

inline void *aligned_malloc(size_t n) {
	return _mm_malloc(n, 64);
}

void linspace(ftype* a, const ftype start, const ftype end, const int n) {
	ftype step = (end - start) / (n - 1);
	ftype value = start;
	for (int i = 0; i < n; ++i) {
		a[i] = value;
		value += step;
	}
}

void dump(ftype* a, int n, const char* s) {
	dprintf("%s", s);
	for (int i = 0; i < n; ++i) {
		dprintf("%.8lf\n", a[i]);
	}
	dprintf("\n");

}

void dump(int* a, int n, const char* s) {
	dprintf("\n%s\n", s);
	for (int i = 0; i < n; ++i) {
		dprintf("%d\n", a[i]);
	}
	dprintf("\n");

}

ftype mean(const ftype data[], const int n) {
	ftype m = 0;
	for (int i = 0; i < n; ++i) {
		m += data[i];
	}
	return m / n;
}

ftype standard_deviation(const ftype data[], const int n, const ftype mean) {
	ftype sum_deviation = 0.0;
	int i;
	for (i = 0; i < n; ++i)
		sum_deviation += (data[i] - mean) * (data[i] - mean);
	return sqrt(sum_deviation / n);
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
