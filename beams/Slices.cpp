/*
 * Slices.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "Slices.h"
#include <algorithm>
#include <math_functions.h>
#include <iterator>
#include <omp.h>

//#include <gsl/gsl_multifit_nlin.h>
//#include <gsl/gsl_vector.h>


Slices::Slices(int _n_slices, int _n_sigma, ftype _cut_left, ftype _cut_right,
               cuts_unit_type _cuts_unit, fit_type _fit_option, bool direct_slicing)
{

   this->n_slices = _n_slices;
   this->cut_left = _cut_left;
   this->cut_right = _cut_right;
   this->cuts_unit = _cuts_unit;
   this->fit_option = _fit_option;
   this->n_sigma = _n_sigma;
   //this->beam_spectrum = 0;
   //this->beam_spectrum_freq = 0;

   //this->h = new ftype[omp_get_num_threads()][n_slices];
   //this->h = (ftype *) malloc(n_threads * n_slices * sizeof(ftype));

   this->n_macroparticles = new ftype[n_slices];
   for (int i = 0; i < n_slices; ++i)
      n_macroparticles[i] = 0;

   this->edges = new ftype[n_slices + 1];
   for (int i = 0; i < n_slices + 1; ++i)
      edges[i] = 0;

   this->bin_centers = new ftype[n_slices];
   for (int i = 0; i < n_slices; ++i)
      bin_centers[i] = 0;

   set_cuts();

   if (direct_slicing)
      track();
}



Slices::~Slices()
{
   util::delete_array(n_macroparticles);
   util::delete_array(bin_centers);
   util::delete_array(edges);
   //delete_array (h);
   //free(h);
}

void Slices::set_cuts()
{
   /*
    *Method to set the self.cut_left and self.cut_right properties. This is
    done as a pre-processing if the mode is set to 'const_space', for
    'const_charge' this is calculated each turn.*

    *The frame is defined by :math:`n\sigma_{RMS}` or manually by the user.
    If not, a default frame consisting of taking the whole bunch +5% of the
    maximum distance between two particles in the bunch will be taken
    in each side of the frame.*
    */
   if (cut_left == 0 && cut_right == 0) {
      if (n_sigma == 0) {
         sort_particles();
         cut_left = Beam->dt[0]
                    - 0.05
                    * (Beam->dt[Beam->n_macroparticles - 1]
                       - Beam->dt[0]);
         cut_right = Beam->dt[Beam->n_macroparticles - 1]
                     + 0.05
                     * (Beam->dt[Beam->n_macroparticles - 1]
                        - Beam->dt[0]);

         //dprintf("cut_left = %e\n", cut_left);
         //dprintf("cut_right = %e\n", cut_right);

      } else {
         ftype mean_coords = mymath::mean(Beam->dt.data(), Beam->n_macroparticles);
         //dprintf("mean coors = %e\n", mean_coords);
         ftype sigma_coords = mymath::standard_deviation(Beam->dt.data(),
                              Beam->n_macroparticles, mean_coords);
         //dprintf("mean coors = %e\n", mean_coords);

         cut_left = mean_coords - n_sigma * sigma_coords / 2;
         cut_right = mean_coords + n_sigma * sigma_coords / 2;

      }
   } else {
      cut_left = convert_coordinates(cut_left, cuts_unit);
      cut_right = convert_coordinates(cut_right, cuts_unit);
   }
   //dprintf("cut_left = %e\n", cut_left);
   //dprintf("cut_right = %e\n", cut_right);

   mymath::linspace(edges, cut_left, cut_right, n_slices + 1);
   //dump(edges, n_slices + 1, "edges\n");
   for (int i = 0; i < n_slices; ++i) {
      bin_centers[i] = (edges[i + 1] + edges[i]) / 2;
   }

}

// TODO not implemented the best way
// If dt, dE and id were in the same struct it would be better
void Slices::sort_particles()
{
   /*
    *Sort the particles with respect to their position.*
    */

   std::sort(&Beam->dE[0], &Beam->dE[Beam->n_macroparticles],
             util::MyComparator(Beam->dt.data()));

   std::sort(&Beam->id[0], &Beam->id[Beam->n_macroparticles],
             util::MyComparator(Beam->dt.data()));
   std::sort(&Beam->dt[0], &Beam->dt[Beam->n_macroparticles],
             util::MyComparator(Beam->dt.data()));

}

inline ftype Slices::convert_coordinates(const ftype cut,
      const cuts_unit_type type)
{
   /*
    *Method to convert a value from one input_unit_type to 's'.*
    */
   if (type == s) {
      return cut;
   } else if (type == rad) {
      return cut / RfP->omega_RF[RfP->counter];
   } else {
      dprintf("WARNING: We were supposed to have either s or rad\n");
      return 0;
   }

}

void Slices::track()
{
   slice_constant_space_histogram();
   if (fit_option == fit_type::gaussian_fit)
      gaussian_fit();
}

/*
void Slices::track(const int start, const int end)
{

   // Track method in order to update the slicing along with the tracker.
   // This will update the beam properties (bunch length obtained from the
   // fit, etc.).

   slice_constant_space_histogram(start, end);
   if (fit_option == fit_type::gaussian_fit)
      gaussian_fit();
}
*/


inline void Slices::slice_constant_space_histogram()
{
   /*
    *Constant space slicing with the built-in numpy histogram function,
    with a constant frame. This gives the same profile as the
    slice_constant_space method, but no compute statistics possibilities
    (the index of the particles is needed).*

    *This method is faster than the classic slice_constant_space method
    for high number of particles (~1e6).*
    */


   histogram(Beam->dt.data(), n_macroparticles, cut_left, cut_right, n_slices,
             Beam->n_macroparticles);

}

/*
inline void Slices::slice_constant_space_histogram(const int start,
      const int end)
{

   // Constant space slicing with the built-in numpy histogram function,
   // with a constant frame. This gives the same profile as the
   // slice_constant_space method, but no compute statistics possibilities
   // (the index of the particles is needed).*

   // This method is faster than the classic slice_constant_space method
   // for high number of particles (~1e6).*


   // TODO why using len and not n_macroparticles?
   // WARNING
   // It is because we remove particles? Then n_macroparticles alive should be used
   // Maybe I need to find a way to re arrange particles
   //int n_threads = omp_get_num_threads();
   int id = omp_get_thread_num();
   //printf("ok here\n");
   //dump(Beam->dt, 10, "beam->dt\n");
   histogram(Beam->dt, &h[id * n_slices], cut_left, cut_right, n_slices,
             Beam->n_macroparticles, start, end);
   //printf("%lf\n", h[id][0]);
   //printf("ok here\n");

   #pragma omp barrier

   int d = 0;
   while ((1 << d) < n_threads) {
      int other = id + (1 << d);
      if (id % (1 << (d + 1)) == 0 && other < n_threads) {
         for (int i = 0; i < n_slices; ++i)
            h[id * n_slices + i] += h[other * n_slices + i];
      }
      d++;
      #pragma omp barrier
   }
   //printf("ok here\n");

   #pragma omp master
   for (int i = 0; i < n_slices; ++i) {
      n_macroparticles[i] = h[i];
   }
   //printf("ok here\n");


}
*/

/*
inline void Slices::histogram(const ftype *__restrict__ input,
                              ftype *__restrict__ output, const ftype cut_left,
                              const ftype cut_right, const int n_slices, const int n_macroparticles,
                              const int start, const int end)
{

//int i;
   //printf("ok here\n");

   ftype a;
   ftype fbin;
   int ffbin;
   const ftype inv_bin_width = n_slices / (cut_right - cut_left);
   //dprintf("inv_bin_width = %e\n", inv_bin_width);
   for (int i = 0; i < n_slices; i++)
      output[i] = 0.0;
   //printf("ok here\n");

   for (int i = start; i < end; i++) {
      a = input[i];
      //printf("a = %e\n", a);
      if ((a < cut_left) || (a > cut_right))
         continue;
      //dprintf("I am here\n");
      fbin = (a - cut_left) * inv_bin_width;
      ffbin = (int)(fbin);
//#pragma omp atomic
      output[ffbin] = output[ffbin] + 1.0;
   }
   //printf("ok here\n");

}
*/

inline void Slices::histogram(const ftype *__restrict__ input,
                              ftype *__restrict__ output, const ftype cut_left,
                              const ftype cut_right, const int n_slices,
                              const int n_macroparticles)
{

   const ftype inv_bin_width = n_slices / (cut_right - cut_left);

   // histogram is faster with ints
   typedef int hist_t;

   //hist_t *res = (hist_t *) calloc(n_slices, sizeof(hist_t));
   //ftype *h = (ftype *) calloc(omp_get_max_threads() * n_slices, sizeof(ftype));
   hist_t *h;
   #pragma omp parallel
   {
      const int threads = omp_get_num_threads();


      const int id = omp_get_thread_num();
      int tile = static_cast<int>((n_macroparticles + threads - 1) / threads);
      int start = id * tile;

      int end = std::min(start + tile, n_macroparticles);
      const int row = id * n_slices;

      #pragma omp single
      h = (hist_t *) calloc(threads * n_slices, sizeof(hist_t));

      for (int i = start; i < end; i++) {
         ftype a = input[i];
         if ((a < cut_left) || (a > cut_right))
            continue;
         int ffbin = static_cast<int>((a - cut_left) * inv_bin_width);
         //h[row + ffbin] = h[row + ffbin] + 1.0;
         h[row + ffbin] = h[row + ffbin] + 1;
      }
      #pragma omp barrier

      tile = (n_slices + threads - 1) / threads;
      start = id * tile;
      end = std::min(start + tile, n_slices);

      for (int i = start; i < end; i++)
         output[i] = 0.0;
      //memset(&output[start], 0, (end-start) * sizeof(ftype));

      for (int i = 0; i < threads; ++i) {
         const int r = i * n_slices;
         for (int j = start; j < end; ++j) {
            //res += h[r + j];
            output[j] += h[r + j];
         }
      }
   }
}

void Slices::track_cuts()
{
   /*
    *Track the slice frame (limits and slice position) as the mean of the
    bunch moves.
    Requires Beam statistics!
    Method to be refined!*
    */
   ftype delta = Beam->mean_dt - 0.5 * (cut_left + cut_right);
   cut_left += delta;
   cut_right += delta;
   for (int i = 0; i < n_slices + 1; ++i) {
      edges[i] += delta;
   }
   for (int i = 0; i < n_slices; ++i) {
      bin_centers[i] += delta;
   }

}

inline void Slices::smooth_histogram(const ftype *__restrict__ input,
                                     ftype *__restrict__ output, const ftype cut_left,
                                     const ftype cut_right, const int n_slices, const int n_macroparticles)
{

   int i;
   ftype a;
   ftype fbin;
   ftype ratioffbin;
   ftype ratiofffbin;
   ftype distToCenter;
   int ffbin = 0;
   int fffbin = 0;
   const ftype inv_bin_width = n_slices / (cut_right - cut_left);
   const ftype bin_width = (cut_right - cut_left) / n_slices;

   for (i = 0; i < n_slices; i++) {
      output[i] = 0.0;
   }

   for (i = 0; i < n_macroparticles; i++) {
      a = input[i];
      if ((a < (cut_left + bin_width * 0.5))
            || (a > (cut_right - bin_width * 0.5)))
         continue;
      fbin = (a - cut_left) * inv_bin_width;
      ffbin = (int)(fbin);
      distToCenter = fbin - (ftype)(ffbin);
      if (distToCenter > 0.5)
         fffbin = (int)(fbin + 1.0);
      ratioffbin = 1.5 - distToCenter;
      ratiofffbin = 1 - ratioffbin;
      if (distToCenter < 0.5)
         fffbin = (int)(fbin - 1.0);
      ratioffbin = 0.5 - distToCenter;
      ratiofffbin = 1 - ratioffbin;
      output[ffbin] = output[ffbin] + ratioffbin;
      output[fffbin] = output[fffbin] + ratiofffbin;
   }
}

void Slices::slice_constant_space_histogram_smooth()
{
   /*
    At the moment 4x slower than slice_constant_space_histogram but smoother.
    */
   smooth_histogram(Beam->dt.data(), this->n_macroparticles, cut_left, cut_right,
                    n_slices, Beam->n_macroparticles);

}

void Slices::rms()
{

   /*
    * Computation of the RMS bunch length and position from the line density
    (bunch length = 4sigma).*
    */
   ftype *lineDenNormalized = new ftype[n_slices];
   ftype *array = new ftype[n_slices];

   ftype timeResolution = bin_centers[1] - bin_centers[0];
   ftype trap = mymath::trapezoid(n_macroparticles, timeResolution,
                                  Beam->n_macroparticles);

   for (int i = 0; i < n_slices; ++i) {
      lineDenNormalized[i] = n_macroparticles[i] / trap;
   }

   for (int i = 0; i < n_slices; ++i) {
      array[i] = bin_centers[i] * lineDenNormalized[i];
   }

   bp_rms = mymath::trapezoid(array, timeResolution, n_slices);

   for (int i = 0; i < n_slices; ++i) {
      array[i] = (bin_centers[i] - bp_rms) * (bin_centers[i] - bp_rms)
                 * lineDenNormalized[i];
   }
   ftype temp = mymath::trapezoid(array, timeResolution, n_slices);
   bl_rms = 4 * sqrt(temp);

   delete[] lineDenNormalized;
   delete[] array;

}

void Slices::fwhm(const ftype shift)
{

   /*
    * Computation of the bunch length and position from the FWHM
    assuming Gaussian line density.*
    */
   int max_i = mymath::max(n_macroparticles, n_slices, 1);
   ftype half_max = shift + 0.5 * (n_macroparticles[max_i] - shift);
   ftype timeResolution = bin_centers[1] - bin_centers[0];
//printf("n_macroparticles.max = %.0lf\n", n_macroparticles[max_i]);
//printf("timeResolution = %e\n", timeResolution);
//printf("half_max = %e\n", half_max);
// First aproximation for the half maximum values

   int i = 0;
   while (n_macroparticles[i] < half_max && i < n_slices)
      i++;
   int taux1 = i;
   i = n_slices - 1;
   while (n_macroparticles[i] < half_max)
      i--;
   int taux2 = i;

//dprintf("taux1, taux2 = %d, %d\n", taux1, taux2);
   ftype t1, t2;

// maybe we could specify what kind of exceptions may occur here
// numerical (divide by zero) or index out of bounds

// TODO something weird is happening here
// Python throws an exception only if you access an element after the end of the array
// but not if you access element before the start of an array
// (in that case it takes the last element of the array)
// Cpp does not throw an exception on eiter occassion
// The right condition is the following in comments
   if (taux1 > 0 && taux2 < n_slices - 1) {
      //if (taux2 < n_slices - 1) {
      try {
         t1 = bin_centers[taux1]
              - (n_macroparticles[taux1] - half_max)
              / (n_macroparticles[taux1]
                 - n_macroparticles[taux1 - 1])
              * timeResolution;
         t2 = bin_centers[taux2]
              + (n_macroparticles[taux2] - half_max)
              / (n_macroparticles[taux2]
                 - n_macroparticles[taux2 + 1])
              * timeResolution;
         //dprintf("t1 = %e\n", t1);
         //dprintf("t2 = %e\n", t2);

         bl_fwhm = 4 * (t2 - t1) / cfwhm;
         bp_fwhm = (t1 + t2) / 2;
         //dprintf("t1, t2 = %e, %e\n", t1, t2);
      } catch (...) {
         bl_fwhm = nan("");
         bp_fwhm = nan("");
      }
   } else {
      //catch (...) {
      //dprintf("taux1, taux2 = %d, %d\n", taux1, taux2);
      //dprintf("exception\n");
      bl_fwhm = nan("");
      bp_fwhm = nan("");
   }

}

ftype Slices::fast_fwhm()
{

   /*
    * Computation of the bunch length and position from the FWHM
    assuming Gaussian line density.*

    height = np.max(self.n_macroparticles)
    index = np.where(self.n_macroparticles > height/2.)[0]
    return cfwhm*(Slices.bin_centers[index[-1]] - Slices.bin_centers[index[0]])
    */
   int max_i = mymath::max(n_macroparticles, Beam->n_macroparticles, 1);
   ftype half_max = 0.5 * n_macroparticles[max_i];

// First aproximation for the half maximum values
// TODO is this correct?

   int i = 0;
   while (n_macroparticles[i] < half_max && i < n_slices)
      i++;
   int taux1 = i;
   i = n_slices - 1;
   while (n_macroparticles[i] < half_max && i >= 0)
      i--;
   int taux2 = i;
   // update bp
   return cfwhm * (bin_centers[taux2] - bin_centers[taux1]);

}

void Slices::fwhm_multibunch()
{
}

void Slices::beam_spectrum_generation(uint n, bool onlyRFFT)
{

   fBeamSpectrumFreq = mymath::rfftfreq(n, bin_centers[1] - bin_centers[0]);

   if (onlyRFFT == false) {
      // TODO remove this when you have moved to vectors
      f_vector_t v(n_macroparticles, n_macroparticles + n_slices);
      //std:: cout << "n is " << n << "\n";
      //std:: cout << "n_slices is " << n_slices << "\n";
      mymath::rfft(v, n, fBeamSpectrum);
      //std:: cout << "size of fBeamSpectrum is " << fBeamSpectrum.size() << "\n";

   }
}

void Slices::beam_profile_derivative()
{
}

void Slices::beam_profile_filter_chebyshev()
{
}

ftype Slices::gauss(const ftype x, const ftype x0, const ftype sx,
                    const ftype A)
{
   return A * exp(-(x - x0) * (x - x0) / 2.0 / (sx * sx));
}

void Slices::gaussian_fit()
{
}

/*

 struct data {
 size_t n;
 double * y;
 };

 int gauss(const gsl_vector * x, void *data, gsl_vector * f) {
 size_t n = ((struct data *) data)->n;
 double *y = ((struct data *) data)->y;

 double A = gsl_vector_get(x, 2);
 double x0 = gsl_vector_get(x, 0);
 double sx = gsl_vector_get(x, 1);

 size_t i;

 for (i = 0; i < n; i++) {
 // Model Yi = A * exp(-lambda * i) + b
 double t = i;
 double Yi = A * exp(-(y[i] - x0) * (y[i] - x0) / 2.0 / (sx * sx));
 gsl_vector_set(f, i, Yi - y[i]);
 }

 return GSL_SUCCESS;
 }

 void Slices::gaussian_fit() {

 // *Gaussian fit of the profile, in order to get the bunch length and
 // position. Returns fit values in units of s.*

 ftype x0, sx, A;
 int max_i = mymath::max(n_macroparticles, n_slices, 1);
 A = n_macroparticles[max_i];
 if (bl_gauss == 0 && bp_gauss == 0) {
 x0 = mymath::mean(beam->dt, beam->n_macroparticles);
 sx = mymath::standard_deviation(beam->dt, beam->n_macroparticles, x0);
 } else {
 x0 = bp_gauss;
 sx = bl_gauss / 4;
 }
 const gsl_multifit_fsolver_type *T; //= gsl_multifit_fdfsolver_lmsder;
 gsl_multifit_fsolver *s;
 int status, info;
 size_t i;
 const size_t n = beam->n_macroparticles;
 const size_t p = 3;

 s = gsl_multifit_fsolver_alloc(T, n, p);
 //gsl_matrix *J = gsl_matrix_alloc(n, p);
 //gsl_matrix *covar = gsl_matrix_alloc(p, p);
 ftype y[n], weights[n];
 struct data d = { n, y };
 gsl_multifit_function f;
 ftype x_init[3] = { x0, sx, A };
 gsl_vector_view x = gsl_vector_view_array(x_init, p);
 gsl_multifit_fsolver_set(s, &f, &x.vector);
 //  TODO experimental values

 status = gsl_multifit_fsolver_driver(s, 20, 100, 1);

 #define FIT(i) gsl_vector_get(s->x, i)
 ftype x0_new = FIT(0);
 ftype sx_new = FIT(1);
 ftype A_new = FIT(2);
 // TODO use gsl to calculate the gaussian fit


 }

 */
