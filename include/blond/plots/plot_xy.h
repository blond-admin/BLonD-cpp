/*
 * plot_xy.h
 *
 *  Created on: Dec 20, 2016
 *      Author: alasheen
 */

#ifndef INCLUDE_BLOND_PLOTS_PLOT_XY_H_
#define INCLUDE_BLOND_PLOTS_PLOT_XY_H_

#include <blond/utilities.h>
namespace blond {
    int plot_xy(f_vector_t xArray, f_vector_t yArray, std::string figureName,
                std::string dirName, std::string saveFigName, double xMax = 0);
} // blond


#endif /* INCLUDE_BLOND_PLOTS_PLOT_XY_H_ */
