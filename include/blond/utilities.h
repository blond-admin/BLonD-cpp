/*
 * utilities.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_UTILITIES_H_
#define INCLUDES_UTILITIES_H_

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <blond/configuration.h>

#ifndef WIN32
#include <mm_malloc.h>
#endif

// DLL export
#if defined(_GCC)
#define EXPORT_DECL __attribute__((visibility("default")))
#define IMPORT_DECL __attribute__((visibility("hidden")))
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32) ||               \
    defined(MSC_VER)
#define EXPORT_DECL __declspec(dllexport)
#define IMPORT_DECL __declspec(dllimport)
#else
#define EXPORT_DECL
#define IMPORT_DECL
#endif

#ifdef SHARED
#define API EXPORT_DECL
#elif STATIC
#define API
#else
#define API IMPORT_DECL
#endif

#include <blond/optionparser.h>

#define dprintf(...) fprintf(stdout, __VA_ARGS__) // Debug printf



#define ALL(c) (c).begin(),(c).end() 
#define FOR(c, it) for(auto it = c.begin(); it != c.end(); it++)

namespace util {

    static inline f_vector_t string_to_double_vector(std::string &numbers)
    {
        f_vector_t res;
        std::istringstream iss(numbers);
        std::copy(std::istream_iterator<double>(iss),
                  std::istream_iterator<double>(),
                  std::back_inserter(res));
        return res;
    }

    static inline std::string double_vector_to_string(f_vector_t &numbers)
    {
        std::ostringstream oss;
        for (auto &n : numbers)
            oss << n << " ";
        return oss.str();
    }

    static inline size_t getFileSize(const std::string &file)
    {
        struct stat st;
        if (stat(file.c_str(), &st) != 0) {
            return 0;
        }
        return st.st_size;
    }

    template <typename T>
    static inline void read_vector_from_file(std::vector<T> &v,
            std::string file)
    {
        v.clear();
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

    template <typename T>
    static inline void dump_to_file(std::vector<T> &v, std::string file)
    {

        std::ofstream output_file(file);

        std::ostream_iterator<T> output_iterator(output_file, "\n");
        std::copy(v.begin(), v.end(), output_iterator);

        output_file.close();
    }

    /*
    template<typename T>
    static inline void dump_to_file(std::vector<std::vector<T>> &v,
                                 std::string file)
    {

     std::ofstream output_file(file);

     std::ostream_iterator<std::string> output_iterator(output_file, "\n");
     for (const auto &v_line : v) {

         for ()
         output_iterator << line;
     }
     // std::copy(v.begin(), v.end(), output_iterator);

     output_file.close();
    }
    */

    static inline char const *GETENV(char const *envstr)
    {
        char const *env = getenv(envstr);
        if (!env)
            return "0";
        else
            return env;
    }

    static inline void *aligned_malloc(size_t n) { return _mm_malloc(n, 64); }

    template <typename T> static inline void delete_array(T *p)
    {
        if (p != NULL)
            delete[] p;
    }

    template <typename T> static inline void zero(T *p, int n)
    {
        for (int i = 0; i < n; ++i) {
            p[i] = 0;
        }
    }

    template <typename T>
    static inline void dump(const T *a, const unsigned n, const char *s)
    {
#ifdef PRINT_RESULTS
        std::cout.precision(PRECISION);
        std::cout << s;
        std::cout << std::scientific << std::showpos;
        for (uint i = 0; i < n; ++i)
            std::cout << a[i] << std::endl;
        std::cout << std::endl;
#endif
    }

    template <typename T>
    static inline void dump(const std::vector<T> &a, const char *s,
                            uint n = 0)
    {
#ifdef PRINT_RESULTS
        n = (n == 0) ? a.size() : n;
        std::cout.precision(PRECISION);
        std::cout << s;
        std::cout << std::scientific << std::showpos;
        for (uint i = 0; i < n; ++i)
            std::cout << a[i] << std::endl;
        std::cout << std::endl;
#endif
    }

    template <typename T> static inline void dump(const T a, const char *s)
    {
#ifdef PRINT_RESULTS
        std::cout.precision(PRECISION);
        std::cout << s;
        std::cout << std::scientific << std::showpos;
        std::cout << a << std::endl;
#endif
    }

    static inline double time_diff(timespec const &end, timespec const &begin)
    {
#ifdef TIMING
        double result;

        result = end.tv_sec - begin.tv_sec;
        result += (end.tv_nsec - begin.tv_nsec) / (double)1000000000;

        return result;
#else
        return 0;
#endif
    }

    static inline void get_time(timespec &ts)
    {
#ifdef TIMING
        auto time = std::chrono::system_clock::now();
        ts.tv_sec = std::chrono::duration_cast<std::chrono::seconds>(
                        time.time_since_epoch())
                    .count();
        ts.tv_nsec = std::chrono::duration_cast<std::chrono::microseconds>(
                         time.time_since_epoch())
                     .count();
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
#ifdef WIN32
        std::shared_ptr<FILE> pipe(_popen(cmd, "r"), _pclose);
#else
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
#endif
        if (!pipe)
            return "ERROR";
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

    struct API Arg : public option::Arg {
        static void printError(const char *msg1, const option::Option &opt,
                               const char *msg2)
        {
            fprintf(stderr, "%s", msg1);
            fwrite(opt.name, opt.namelen, 1, stderr);
            fprintf(stderr, "%s", msg2);
        }

        static option::ArgStatus Unknown(const option::Option &option,
                                         bool msg)
        {
            if (msg)
                printError("Unknown option '", option, "'\n");
            return option::ARG_ILLEGAL;
        }

        static option::ArgStatus Required(const option::Option &option,
                                          bool msg)
        {
            if (option.arg != 0)
                return option::ARG_OK;

            if (msg)
                printError("Option '", option, "' requires an argument\n");
            return option::ARG_ILLEGAL;
        }

        static option::ArgStatus NonEmpty(const option::Option &option,
                                          bool msg)
        {
            if (option.arg != 0 && option.arg[0] != 0)
                return option::ARG_OK;

            if (msg)
                printError("Option '", option,
                           "' requires a non-empty argument\n");
            return option::ARG_ILLEGAL;
        }

        static option::ArgStatus Numeric(const option::Option &option,
                                         bool msg)
        {
            // printf("Inside here\n");
            char *endptr = 0;
            if (option.arg != 0 && strtol(option.arg, &endptr, 10)) {
            };
            if (endptr != option.arg && *endptr == 0)
                return option::ARG_OK;

            if (msg)
                printError("Option '", option,
                           "' requires a numeric argument\n");
            return option::ARG_ILLEGAL;
        }
    };
}

#endif /* INCLUDES_UTILITIES_H_ */
