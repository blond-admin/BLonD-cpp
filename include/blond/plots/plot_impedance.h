

#ifndef INCLUDE_BLOND_PLOTS_PLOT_IMPEDANCE_H_
#define INCLUDE_BLOND_PLOTS_PLOT_IMPEDANCE_H_

#include <blond/impedances/InducedVoltage.h>
#include <blond/beams/Slices.h>
#include <blond/python.h>
namespace blond {
int plot_impedance_vs_frequency(int counter, InducedVoltageFreq *indVoltFreq,
                                 Slices *slices,
                                 std::string option1 = "sum",
                                 std::string option2 = "no_spectrum",
                                 std::string option3 = "freq_fft",
                                 std::string style = "-",
                                 f_vector_t cut_left_right = {},
                                 f_vector_t cut_up_down = {},
                                 std::string dirname = "fig");



int plot_induced_voltage_vs_bin_centers(int counter,
        TotalInducedVoltage *totIndVolt,
        Slices *slices,
        std::string style = "-",
        std::string dirname = "fig");
} // blond

#endif /* INCLUDE_BLOND_PLOTS_PLOT_IMPEDANCE_H_ */
