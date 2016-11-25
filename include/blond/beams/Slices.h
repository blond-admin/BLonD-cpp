/*
 * Slices.h
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_SLICES_H_
#define BEAMS_SLICES_H_

class Slices;

#include <blond/configuration.h>
#include <blond/utilities.h>
#include <blond/input_parameters/RfParameters.h>
#include <map>


class API Slices {
private:
    const double cfwhm = 2 * std::sqrt(2 * std::log(2));

    f_vector_t gaussian_filter1d(f_vector_t &x, int sigma,
                                 int order, std::string mode);
    f_vector_t gradient(f_vector_t &x, double dist);
public:
    enum cuts_unit_t { s, rad };
    enum fit_t { normal, gaussian };

    Beams *beam;
    RfParameters *rfp;

    double bl_fwhm, bp_fwhm;
    double bp_rms, bl_rms;
    int n_slices;
    double cut_left;
    double cut_right;
    int n_sigma;
    cuts_unit_t cuts_unit;
    f_vector_t n_macroparticles;
    f_vector_t edges;
    f_vector_t bin_centers;
    fit_t fit_option;
    complex_vector_t fBeamSpectrum;
    f_vector_t fBeamSpectrumFreq;
    double bl_gauss;
    double bp_gauss;

    Slices(RfParameters *RfP, Beams *Beam,
           int _n_slices, int _n_sigma = 0, double cut_left = 0,
           double cut_right = 0, cuts_unit_t cuts_unit = s,
           fit_t fit_option = normal, bool direct_slicing = false);

    ~Slices();
    void track();
    // double fast_fwhm();
    void fwhm(const double shift = 0);
    void beam_spectrum_generation(int n, bool onlyRFFT = false);
    void beam_profile_derivative(f_vector_t &x,
                                 f_vector_t &derivative,
                                 std::string mode = "gradient");

    void beam_profile_filter_chebyshev(std::map<std::string, std::string>
                                       filter_option,
                                       int &nCoefficients,
                                       f_vector_t &b,
                                       f_vector_t &a);

    void beam_profile_filter_chebyshev(std::map<std::string, std::string>
                                       filter_option,
                                       int &nCoefficients,
                                       f_vector_t &transferFreq,
                                       complex_vector_t &transferGain);
    void set_cuts();
    void sort_particles();
    double convert_coordinates(double cut, cuts_unit_t type);

    void histogram(const double *__restrict input, double *__restrict output,
                   const double cut_left, const double cut_right,
                   const int n_slices, const int n_macroparticles);

    void smooth_histogram(const double *__restrict input,
                          double *__restrict output, const double cut_left,
                          const double cut_right, const int n_slices,
                          const int n_macroparticles);
    void slice_constant_space_histogram();
    void track_cuts();
    void slice_constant_space_histogram_smooth();
    void rms();
    void gaussian_fit();

    void fwhm_multibunch();
    // when intensity effects
};
#endif /* BEAMS_SLICES_H_ */
