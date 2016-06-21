/*
 * Slices.h
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_SLICES_H_
#define BEAMS_SLICES_H_

class Slices;

#include <blond/beams/Beams.h>
#include <blond/configuration.h>
#include <blond/constants.h>
#include <blond/fft.h>
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/utilities.h>

const ftype cfwhm = 2 * sqrt(2 * log(2));

enum cuts_unit_type { s, rad };

enum fit_type { normal_fit, gaussian_fit };

class Slices {
  public:
    Slices(int _n_slices, int _n_sigma = 0, ftype cut_left = 0,
           ftype cut_right = 0, cuts_unit_type cuts_unit = s,
           fit_type fit_option = normal_fit, bool direct_slicing = false);

    ~Slices();

    // void track(const int start, const int end);
    void track();

    // void zero_histogram();
    ftype fast_fwhm();

    void fwhm(const ftype shift = 0);

    ftype bl_fwhm, bp_fwhm;
    ftype bp_rms, bl_rms;
    int n_slices;
    ftype cut_left;
    ftype cut_right;
    int n_sigma;
    cuts_unit_type cuts_unit;
    ftype* n_macroparticles;
    ftype* edges;
    ftype* bin_centers;
    fit_type fit_option;

    complex_vector_t fBeamSpectrum;

    f_vector_t fBeamSpectrumFreq;

    ftype bl_gauss = 0;
    ftype bp_gauss = 0;

    void beam_spectrum_generation(uint n, bool onlyRFFT = false);

    void beam_profile_derivative();

    void beam_profile_filter_chebyshev();

  private:
    // ftype *h;
    void set_cuts();

    void sort_particles();

    inline ftype convert_coordinates(ftype cut, cuts_unit_type type);

    // inline void histogram(const ftype *__restrict__ input,
    //                      ftype *__restrict__ output, const ftype cut_left,
    //                      const ftype cut_right, const int n_slices,
    //                      const int n_macroparticles, const int start, const
    //                      int end);
    inline void histogram(const ftype* __restrict__ input,
                          ftype* __restrict__ output, const ftype cut_left,
                          const ftype cut_right, const int n_slices,
                          const int n_macroparticles);

    inline void smooth_histogram(const ftype* __restrict__ input,
                                 ftype* __restrict__ output,
                                 const ftype cut_left, const ftype cut_right,
                                 const int n_slices,
                                 const int n_macroparticles);

    // inline void slice_constant_space_histogram(const int start, const int
    // end);
    inline void slice_constant_space_histogram();

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
