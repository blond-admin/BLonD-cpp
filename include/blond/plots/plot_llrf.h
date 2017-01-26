

#ifndef INCLUDE_BLOND_PLOTS_PLOT_LLRF_H_
#define INCLUDE_BLOND_PLOTS_PLOT_LLRF_H_

#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/python.h>


namespace blond {
    int plot_noise_spectrum(f_vector_t &frequency, f_vector_t &spectrum,
                            int sampling = 1, std::string dirname = "fig",
                            int figno = 0);


    int plot_phase_noise(f_vector_t &frequency, f_vector_t &spectrum,
                         int sampling = 1, std::string dirname = "fig",
                         int figno = 0);


    int plot_PL_bunch_phase(RfParameters *RfP, std::string h5data,
                            int output_freq = 1, std::string dirname = "fig");


    int plot_PL_RF_phase(RfParameters *RfP, std::string h5data,
                         int output_freq = 1, std::string dirname = "fig");


    int plot_PL_phase_corr(RfParameters *RfP, std::string h5data,
                           int output_freq = 1, std::string dirname = "fig");


    int plot_PL_RF_freq(RfParameters *RfP, std::string h5data,
                        int output_freq = 1, std::string dirname = "fig");


    int plot_PL_freq_corr(RfParameters *RfP, std::string h5data,
                          int output_freq = 1, std::string dirname = "fig");


    int plot_RF_phase_error(RfParameters *RfP, std::string h5data,
                            int output_freq = 1, std::string dirname = "fig");


    int plot_RL_radial_error(RfParameters *RfP, std::string h5data,
                             int output_freq = 1, std::string dirname = "fig");


    int plot_COM_motion(RfParameters *RfP,
                        std::string h5data, int output_freq = 1,
                        std::string dirname = "fig");


    int plot_LHCNoiseFB(RfParameters *RfP, std::string h5data,
                        int output_freq = 1, std::string dirname = "fig");



    int plot_LHCNoiseFB_FWHM(RfParameters *RfP,
                             std::string h5data, int output_freq = 1,
                             std::string dirname = "fig");



    int plot_LHCNoiseFB_FWHM_bbb(RfParameters *RfP,
                                 std::string h5data, int output_freq = 1,
                                 std::string dirname = "fig");

} // blond
#endif /* INCLUDE_BLOND_PLOTS_PLOT_LLRF_H_ */
