

#ifndef INCLUDE_BLOND_PLOTS_PLOT_SLICES_H_
#define INCLUDE_BLOND_PLOTS_PLOT_SLICES_H_

#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Slices.h>


void plot_beam_profile(Slices *Slices,
                       int counter,
                       std::string style = "-",
                       std::string dirname = "fig");

void plot_beam_profile_derivative(Slices *Slices,
                                  int counter,
                                  std::string style = "-",
                                  std::string dirname = "fig",
                                  int_vector_t numbers = int_vector_t({3}));

void plot_beam_spectrum(Slices *Slices, int counter,
                        std::string style = "-",
                        std::string dirname = "fig");



#endif /* INCLUDE_BLOND_PLOTS_PLOT_SLICES_H_ */
