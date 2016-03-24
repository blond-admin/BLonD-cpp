/*
 * utilities.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "utilities.h"

char const* GETENV(char const* envstr) {
	char const* env = getenv(envstr);
	if (!env)
		return "0";
	else
		return env;
}

void *aligned_malloc(size_t n) {
	return _mm_malloc(n, 64);
}

void delete_array(ftype *p) {
	if (p != NULL)
		delete[] p;
}

void delete_array(bool *p) {
	if (p != NULL)
		delete[] p;
}

void delete_array(int *p) {
	if (p != NULL)
		delete[] p;
}

void dump(ftype* a, int n, const char* s) {
	dprintf("%s", s);
	for (int i = 0; i < n; ++i) {
		dprintf("%.12lf\n", a[i]);
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

double time_diff(timespec const& end, timespec const& begin) {
#ifdef TIMING
	double result;

	result = end.tv_sec - begin.tv_sec;
	result += (end.tv_nsec - begin.tv_nsec) / (double) 1000000000;

	return result;
#else
	return 0;
#endif
}

void get_time(timespec& ts) {
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

timespec get_time() {
	timespec t;
#ifdef TIMING
	get_time(t);
#endif
	return t;
}

double time_elapsed(timespec const& begin) {
#ifdef TIMING
	timespec now;
	get_time(now);
	return time_diff(now, begin);
#else
	return 0;
#endif
}

void print_time(char const* prompt, timespec const& begin,
		timespec const& end) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, time_diff(end, begin));
#endif
}

void print_time(char const* prompt, double diff) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, diff);
#endif
}

void print_time_elapsed(char const* prompt, timespec const& begin) {
#ifdef TIMING
	dprintf("%s : %.3f\n", prompt, time_elapsed(begin));
#endif
}
