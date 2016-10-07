/*
 * Slices.h
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_SLICES_H_
#define BEAMS_SLICES_H_

#include <blond/configuration.h>
#include <blond/utilities.h>

const ftype cfwhm = 2 * sqrt(2 * log(2));



class API Slices {
private:

    f_vector_t gaussian_filter1d(f_vector_t &x, int sigma,
                                 int order, std::string mode);
    f_vector_t gradient(f_vector_t &x, ftype dist);
public:
    enum cuts_unit_t { s, rad };
    enum fit_t { normal, gaussian };

    ftype bl_fwhm, bp_fwhm;
    ftype bp_rms, bl_rms;
    int n_slices;
    ftype cut_left;
    ftype cut_right;
    int n_sigma;
    cuts_unit_t cuts_unit;
    int_vector_t n_macroparticles;
    f_vector_t fFloatMacroparticles;
    f_vector_t edges;
    f_vector_t bin_centers;
    fit_t fit_option;
    complex_vector_t fBeamSpectrum;
    f_vector_t fBeamSpectrumFreq;
    ftype bl_gauss = 0;
    ftype bp_gauss = 0;

    Slices(uint _n_slices, int _n_sigma = 0, ftype cut_left = 0,
           ftype cut_right = 0, cuts_unit_t cuts_unit = s,
           fit_t fit_option = normal, bool direct_slicing = false);

    ~Slices();
    void track();
    // ftype fast_fwhm();
    void fwhm(const ftype shift = 0);
    void beam_spectrum_generation(uint n, bool onlyRFFT = false);
    void beam_profile_derivative(f_vector_t &x,
                                 f_vector_t &derivative,
                                 std::string mode = "gradient");
    void beam_profile_filter_chebyshev();
    void set_cuts();
    void sort_particles();
    ftype convert_coordinates(ftype cut, cuts_unit_t type);

    void histogram(const ftype *__restrict input, int *__restrict output,
                   const ftype cut_left, const ftype cut_right,
                   const int n_slices, const int n_macroparticles);
    void smooth_histogram(const ftype *__restrict input,
                          ftype *__restrict output, const ftype cut_left,
                          const ftype cut_right, const int n_slices,
                          const int n_macroparticles);
    void slice_constant_space_histogram();
    void track_cuts();
    void slice_constant_space_histogram_smooth();
    void rms();
    ftype gauss(const ftype x, const ftype x0, const ftype sx, const ftype A);

    // not for now
    void gaussian_fit();
    // not for now
    void fwhm_multibunch();
    // when intensity effects
};
#endif /* BEAMS_SLICES_H_ */
