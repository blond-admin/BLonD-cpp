


#ifndef INCLUDE_BLOND_PLOTS_PLOT_PARAMETERS_H_
#define INCLUDE_BLOND_PLOTS_PLOT_PARAMETERS_H_

#include <blond/configuration.h>
#include <blond/python.h>
#include <string>
namespace blond {
    int plot_voltage_programme(f_vector_t &time, f_vector_t &voltage,
                               int sampling = 1, std::string dirname = "fig",
                               int figno = 0);

} // blond


#endif /* INCLUDE_BLOND_PLOTS_PLOT_PARAMETERS_H_ */
