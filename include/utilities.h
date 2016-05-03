/*
 * utilities.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_UTILITIES_H_
#define INCLUDES_UTILITIES_H_

#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <mm_malloc.h>
#include <sys/time.h>
#include "configuration.h"
#include <memory>
#include <optionparser.h>

#define dprintf(...)    fprintf(stdout, __VA_ARGS__)     // Debug printf

#define CHECK_ERROR(a)                                       \
   if (a)                                                    \
   {                                                         \
      perror("Error at line\n\t" #a "\nSystem Msg");         \
      assert ((a) == 0);                                     \
   }

namespace util {

   template<typename T>
   static inline void read_vector_from_file(std::vector<T> &v, std::string file)
   {

      std::ifstream source(file);
      if (!source.good()) {
         std::cout << "Error: file " << file << " does not exist\n";
         source.close();
         exit(-1);
      }

      for (std::string line; std::getline(source, line);) {
         std::istringstream in(line);
         T type;
         while (in >> type)
            v.push_back(type);
      }

      source.close();
   }

// sort an array with regards to another array
   struct MyComparator {
      ftype *a;
      MyComparator(ftype *_a) :
         a(_a)
      {
      }

      bool operator()(ftype i1, ftype i2)
      {
         return i1 < i2;
      }
   };

   static inline char const *GETENV(char const *envstr)
   {
      char const *env = getenv(envstr);
      if (!env)
         return "0";
      else
         return env;
   }

   static inline void *aligned_malloc(size_t n)
   {
      return _mm_malloc(n, 64);
   }

   template<typename T>
   static inline void delete_array(T *p)
   {
      if (p != NULL)
         delete[] p;
   }

   template<typename T>
   static inline void zero(T *p, int n)
   {
      for (int i = 0; i < n; ++i) {
         p[i] = 0;
      }
   }

   static inline void dump(ftype *a, int n, const char *s)
   {

#ifdef PRINT_RESULTS
      dprintf("%s", s);
      for (int i = 0; i < n; ++i) {
         dprintf("%+.6e\n", a[i]);
      }
      dprintf("\n");
#endif

   }

   static inline void dump(ftype a, const char *s)
   {

#ifdef PRINT_RESULTS

      dprintf("%s: %+.6e\n", s, a);
#endif

   }


   static inline void dump(int *a, int n, const char *s)
   {

#ifdef PRINT_RESULTS
      dprintf("%s", s);
      for (int i = 0; i < n; ++i) {
         dprintf("%+d\n", a[i]);
      }
      dprintf("\n");
#endif


   }

   static inline void dump(int a, const char *s)
   {

#ifdef PRINT_RESULTS
      dprintf("%s: %+d\n", s, a);
#endif

   }

   static inline double time_diff(timespec const &end, timespec const &begin)
   {
#ifdef TIMING
      double result;

      result = end.tv_sec - begin.tv_sec;
      result += (end.tv_nsec - begin.tv_nsec) / (double) 1000000000;

      return result;
#else
      return 0;
#endif
   }

   static inline void get_time(timespec &ts)
   {

#ifdef TIMING
      struct timeval tv;
      gettimeofday(&tv, NULL);
      ts.tv_sec = tv.tv_sec;
      ts.tv_nsec = tv.tv_usec * 1000;
#endif
   }

   static inline timespec get_time()
   {
      timespec t;
#ifdef TIMING
      get_time(t);
#endif
      return t;
   }

   static inline double time_elapsed(timespec const &begin)
   {
#ifdef TIMING
      timespec now;
      get_time(now);
      return time_diff(now, begin);
#else
      return 0;
#endif
   }

   static inline void print_time(char const *prompt, timespec const &begin,
                                 timespec const &end)
   {
#ifdef TIMING
      dprintf("%s : %.3f\n", prompt, time_diff(end, begin));
#endif
   }

   static inline void print_time(char const *prompt, double diff)
   {
#ifdef TIMING
      dprintf("%s : %.3f\n", prompt, diff);
#endif
   }

   static inline void print_time_elapsed(char const *prompt,
                                         timespec const &begin)
   {
#ifdef TIMING
      dprintf("%s : %.3f\n", prompt, time_elapsed(begin));
#endif
   }


   static inline std::string exec(const char *cmd)
   {
      std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
      if (!pipe) return "ERROR";
      char buffer[128];
      std::string result = "";
      while (!feof(pipe.get())) {
         if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
      }
      return result;
   }

   static inline std::string read_from_file(std::string filename)
   {
      std::ifstream t(filename);
      std::stringstream res;
      res << t.rdbuf();
      return res.str();
   }

   struct Arg: public option::Arg {
      static void printError(const char *msg1, const option::Option &opt, const char *msg2)
      {
         fprintf(stderr, "%s", msg1);
         fwrite(opt.name, opt.namelen, 1, stderr);
         fprintf(stderr, "%s", msg2);
      }

      static option::ArgStatus Unknown(const option::Option &option, bool msg)
      {
         if (msg) printError("Unknown option '", option, "'\n");
         return option::ARG_ILLEGAL;
      }

      static option::ArgStatus Required(const option::Option &option, bool msg)
      {
         if (option.arg != 0)
            return option::ARG_OK;

         if (msg) printError("Option '", option, "' requires an argument\n");
         return option::ARG_ILLEGAL;
      }

      static option::ArgStatus NonEmpty(const option::Option &option, bool msg)
      {
         if (option.arg != 0 && option.arg[0] != 0)
            return option::ARG_OK;

         if (msg) printError("Option '", option, "' requires a non-empty argument\n");
         return option::ARG_ILLEGAL;
      }

      static option::ArgStatus Numeric(const option::Option &option, bool msg)
      {
         //printf("Inside here\n");
         char *endptr = 0;
         if (option.arg != 0 && strtol(option.arg, &endptr, 10)) {};
         if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;

         if (msg) printError("Option '", option, "' requires a numeric argument\n");
         return option::ARG_ILLEGAL;
      }
   };


}

#endif /* INCLUDES_UTILITIES_H_ */
