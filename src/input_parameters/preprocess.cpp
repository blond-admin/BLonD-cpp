#include <blond/input_parameters/preprocess.h>
using namespace blond;

void preprocess_ramp(GeneralParameters::particle_t particle_type,
                     double circumference, f_vector_t &time, f_vector_t &data,
                     f_vector_t &time_interp, f_vector_t &momentum_interp,
                     std::string data_type, std::string interpolation,
                     double smoothing, int flat_bottom, int t_start,
                     int t_end, bool plot, std::string figdir,
                     std::string figname, int sampling, double user_mass)
{}

void preprocess_rf_params(GeneralParameters *GP, f_vector_2d_t &time_arrays,
                          f_vector_2d_t &data_arrays, f_vector_t &data_interp,
                          std::string interpolation,
                          bool plot, std::string figdir,
                          std::vector<std::string> figname, int sampling)
{}

