#ifndef INCLUDE_BLOND_INPUT_PARAMETERS_PREPROCESS_H_
#define INCLUDE_BLOND_INPUT_PARAMETERS_PREPROCESS_H_

#include <blond/GeneralParameters.h>

void preprocess_ramp(GeneralParameters::particle_t particle_type,
                     double circumference, f_vector_t &time, f_vector_t &data,
                     std::string data_type = "momentum",
                     std::string interpolation = "linear",
                     int smoothing = 0, int flat_bottom = 0, int t_start = 0,
                     int t_end = -1, bool plot = false,
                     std::string figdir = "fig", std::string figname = "data",
                     int sampling = 1, double user_mass = 0);
void preprocess_rf_params();

#endif /* INCLUDE_BLOND_INPUT_PARAMETERS_PREPROCESS_H_ */