

#ifndef INCLUDE_BLOND_PLOTS_PLOT_BEAMS_H_
#define INCLUDE_BLOND_PLOTS_PLOT_BEAMS_H_

#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Beams.h>
#include <blond/beams/Slices.h>


void plot_long_phase_space(GeneralParameters *GP, RfParameters *RfP,
                           Beams *Beam, double xmin, double xmax,
                           double ymin, double ymax, std::string xunit = "s",
                           int sampling = 1, bool separatrix_plot = false,
                           bool histograms_plot = true, std::string dirname = "fig",
                           int alpha = 1);


void plot_bunch_length_evol(RfParameters *RfP, std::string h5data,
                            int output_freq = 1, std::string dirname = "fig");


void plot_bunch_length_evol_gaussian(RfParameters *RfP,
                                     std::string h5data, int output_freq = 1,
                                     std::string dirname = "fig");


void plot_position_evol(RfParameters *RfP, std::string h5data,
                        int output_freq = 1, std::string style = ".",
                        std::string dirname = "fig");


void plot_energy_evol(RfParameters *RfP, std::string h5data,
                      int output_freq = 1, std::string style = ".",
                      std::string dirname = "fig");


void plot_transmitted_particles(RfParameters *RfP, std::string h5data,
                                int output_freq = 1, std::string style = ".",
                                std::string dirname = "fig");


#endif /* INCLUDE_BLOND_PLOTS_PLOT_BEAMS_H_ */
