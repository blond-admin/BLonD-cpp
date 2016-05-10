/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 15:25:33
*/

#include "Intensity.h"
#include "utilities.h"
#include "constants.h"


Resonators::Resonators(std::vector<ftype> &RS, std::vector<ftype> &FrequencyR,
                       std::vector<ftype> &Q)
{

   /*
   *Impedance contribution from resonators, analytic formulas for both wake and impedance. The resonant modes (and the corresponding R and Q)
   can be inputed as a list in case of several modes.*

   *The model is the following:*

   .. math::

     Z(f) = \\frac{R}{1 + j Q \\left(\\frac{f}{f_r}-\\frac{f_r}{f}\\right)}

   .. math::

     W(t>0) = 2\\alpha R e^{-\\alpha t}\\left(\\cos{\\bar{\\omega}t} - \\frac{\\alpha}{\\bar{\\omega}}\\sin{\\bar{\\omega}t}\\right)

     W(0) = \\alpha R

   .. math::

     \\omega_r = 2 \\pi f_r

     \\alpha = \\frac{\\omega_r}{2Q}

     \\bar{\\omega} = \\sqrt{\\omega_r^2 - \\alpha^2}
   */


   fRS = RS;
   fFrequencyR = FrequencyR;
   fQ = Q;
   fNResonators = RS.size();

   for (unsigned int i = 0; i < fNResonators; ++i) {
      fOmegaR.push_back(2 * constant::pi * fFrequencyR[i]);
   }

   /*
   fTimeArray = NULL;
   fFreqArray = NULL;
   fWake = NULL;
   fImpedance = NULL;
   */
}



void Resonators::wake_calc(std::vector<ftype> &NewTimeArray)
{
   /*
   * Wake calculation method as a function of time.*
   */
   fTimeArray = NewTimeArray;
   fWake.resize(fTimeArray.size());
   std::fill_n(fWake.begin(), fWake.size(), 0);
   //util::dump(&fWake[0], 10, "start wake ");

   for (uint i = 0; i < fNResonators; ++i) {
      ftype alpha = fOmegaR[i] / (2 * fQ[i]);
      ftype omega_bar = std::sqrt(fOmegaR[i] * fOmegaR[i] -
                                  alpha * alpha);
      //util::dump(&alpha, 1, "alpha ");
      //util::dump(&omega_bar, 1, "omega_bar ");

      for (uint j = 0; j < fWake.size(); ++j) {
         ftype temp = fTimeArray[j];
         //util::dump(&temp, 1, "temp ");
         int sign = (temp > 0) - (temp < 0);
         fWake[j]+= (sign + 1) * fRS[i] * alpha *
                     std::exp(-alpha * temp) *
                     (std::cos(omega_bar * temp) -
                      alpha / omega_bar * std::sin(omega_bar * temp));
      }
   }

}


void Resonators::imped_calc(std::vector<ftype> &NewFrequencyArray)
{
   /*
   * Impedance calculation method as a function of frequency.*
   */

   fFreqArray = NewFrequencyArray;
   fImpedance.resize(fFreqArray.size());
   //fImpedance[0] = complex_t(0, 0);
   std::fill_n(fImpedance.begin(), fImpedance.size(), complex_t(0, 0));
   //std::cout << complex_t(1,0) / complex_t(1,1) << '\n';
   for (uint i = 0; i < fNResonators; ++i) {
      for (uint j = 1; j < fImpedance.size(); ++j) {
         fImpedance[j] += complex_t(fRS[i], 0) /
                          complex_t(1, fQ[i] *
                                    (fFreqArray[j] / fFrequencyR[i] -
                                     fFrequencyR[i] / fFreqArray[j]));
         /*
         if(i==0 && j < 10){
            std::cout << complex_t(fRS[i], 0) << "\n";
            std::cout << complex_t(1, fQ[i] *
                                    (fFreqArray[j] / fFrequencyR[i] -
                                     fFrequencyR[i] / fFreqArray[j])) << "\n";
            std::cout << complex_t(fRS[i], 0) / complex_t(1, fQ[i] *
                                    (fFreqArray[j] / fFrequencyR[i] -
                                     fFrequencyR[i] / fFreqArray[j])) << "\n";
         }
         */
      }
   }


}

InputTable::InputTable()
{

}


void InputTable::wake_calc(std::vector<ftype> &NewTimeArray)
{
}

void InputTable::imped_calc(std::vector<ftype> &NewFrequencyArray)
{
}

