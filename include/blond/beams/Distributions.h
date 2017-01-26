/*
 * Distributions.h
 *
 *  Created on: Mar 16, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_DISTRIBUTIONS_H_
#define BEAMS_DISTRIBUTIONS_H_

#include <blond/configuration.h>
#include <blond/python.h>
#include <blond/beams/Beams.h>
#include <blond/trackers/Tracker.h>
#include <blond/impedances/InducedVoltage.h>
#include <map>
#include <functional>
namespace blond {

    struct line_density_t {
        f_vector_t hamiltonian_coord;
        f_vector_t density_function;
        f_vector_t time_line_den;
        f_vector_t line_density;
        line_density_t(const f_vector_t &ham_c,
                       const f_vector_t &den_func,
                       const f_vector_t &time_l,
                       const f_vector_t &line_den) :
            hamiltonian_coord(ham_c),
            density_function(den_func),
            time_line_den(time_l),
            line_density(line_den)
        {}
    };


    struct distribution_denstity_t {
        f_vector_t time_coord_low_res;
        f_vector_t line_density;
        distribution_denstity_t(const f_vector_t &time_coord,
                                const f_vector_t &line_density) :
            time_coord_low_res(time_coord),
            line_density(line_density)
        {}
    };

    struct multi_t {
        double d;
        std::string s;
        int i;
        f_vector_t v;
        std::function<f_vector_t(const f_vector_t &, std::string,
                                 double, double)> f;
        // std::function<double(double)> f;
        multi_t() {}
        multi_t(double _d) : d(_d) {}
        multi_t(std::string _s) : s(_s) {}
        multi_t(int _i) : i(_i) {}
        multi_t(const f_vector_t &_v) : v(_v) {}
        multi_t(std::function<f_vector_t(const f_vector_t &, std::string,
                                         double, double)> _f) : f(_f) {}
        // multi_t(std::function<double(double)> _f) : f(_f) {}

    };


    line_density_t
    matched_from_line_density(Beams *beam,
                              FullRingAndRf *full_ring,
                              std::map<std::string, std::string> line_density_opt,
                              FullRingAndRf::main_harmonic_t main_harmonic_opt =
                                  FullRingAndRf::lowest_freq,
                              TotalInducedVoltage *totVolt = nullptr,
                              std::string plot = "",
                              std::string figdir = "fig",
                              std::string half_option = "first",
                              std::map<std::string, f_vector_t> extraVoltageDict =
                                  std::map<std::string, f_vector_t>(),
                              int n_iterations_input = 100,
                              int seed = 0);


    distribution_denstity_t
    matched_from_distribution_density(Beams *beam,
                                      FullRingAndRf *full_ring,
                                      std::map<std::string, multi_t> distribution_opt,
                                      FullRingAndRf::main_harmonic_t main_harmonic_opt =
                                          FullRingAndRf::lowest_freq,
                                      TotalInducedVoltage *totVolt = nullptr,
                                      std::map<std::string, f_vector_t> extraVoltageDict =
                                          std::map<std::string, f_vector_t>(),
                                      int n_iterations_input = 1,
                                      int seed = 0);


    /*
    void matched_from_distribution_density(FullRingAndRf *full_ring,
                                           std::map<std::string, std::string> distribution_opt,
                                           std::string main_harmonic = "lowest_freq",
                                           int n_iterations_input = 1,
                                           std::map<std::string, std::string> extraVoltageDict =
                                                   std::map<std::string, std::string>(),
                                           int seed = 0);
    */

    void longitudinal_bigaussian(GeneralParameters *GP, RfParameters *RfP,
                                 Beams *Beam, double sigma_dt, double sigma_dE = 0,
                                 int seed = 0, bool reinsertion = false);


    f_vector_t distribution_density_function(const f_vector_t &action_array,
            const std::string &dist_type, const double length, double exponent = 0.);


    f_vector_2d_t distribution_density_function(const f_vector_2d_t &action_array,
            const std::string &dist_type, const double length, double exponent = 0.0);

    f_vector_t line_density_function(const f_vector_t &coord_array,
                                     const std::string &dist_type,
                                     const double bunch_length,
                                     const double bunch_position = 0.,
                                     double exponent = 0.);

    void plot_generated_bunch(f_vector_t &time_line_den, f_vector_t &line_density,
                              f_vector_t &time_coord_for_grid,
                              f_vector_t &rec_line_den,
                              std::string plot, std::string figdir);


} // blond
#endif /* BEAMS_DISTRIBUTIONS_H_ */
