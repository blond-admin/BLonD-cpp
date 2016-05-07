/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 15:25:33
*/

#include "Intensity.h"
#include "utilities.h"

Resonators::Resonators(uint NResonators, std::vector<ftype> RS,
                       std::vector<ftype> FrequencyR, std::vector<ftype> Q)
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
   fNResonators = NResonators;

   for (int i = 0; i < n_resonators; ++i) {
      fOmegaR.push_back(2 * constants::pi * fFrequencyR[i]);
   }

   /*
   fTimeArray = NULL;
   fFreqArray = NULL;
   fWake = NULL;
   fImpedance = NULL;
   */
}


Resonators::~Resonators()
{
   /*
   util::delete_array(fTimeArray);
   util::delete_array(fFreqArray);
   util::delete_array(fWake);
   util::delete_array(fImpledance);
   util::delete_array(fRS);
   util::delete_array(fFrequencyR);
   util::delete_array(fNResonators);
   util::delete_array(fOmegaR);
   */
}

void Resonators::wake_calc(std::vector<ftype> NewTimeArray)
{
   /*
   * Wake calculation method as a function of time.*
   */
   fTimeArray = NewTimeArray;
   fWake.resize(fTimeArray.size());
   //std::fill_n(fWake, fNResonators, 0);

   for (int i = 0; i < fNResonators; ++i) {
      ftype alpha = fOmegaR[i] / (2 * fQ[i]);
      ftype omega_bar = sqrt(fOmegaR[i] * fOmegaR[i] - alpha * alpha);
      for (int j = 0; j < fWake.size(); ++j) {
         ftype temp = fTimeArray[j];
         fWake[i] += (signbit(temp) + 1) * fRS[i] * alpha *
                     exp(-alpha * temp) * (cos(omega_bar * temp)
                     - alpha / omega_bar * sin(omega_bar * temp));
      }
   }

}


void Resonators::imped_calc(std::vector<ftype> NewFrequencyArray)
{
   /*
   * Impedance calculation method as a function of frequency.*
   */

   fFreqArray = NewFrequencyArray;
   // TODO
   // figure out how to use complex numbers here

}

InputTable::InputTable() {}
InputTable::~InputTable() {}

void InputTable::wake_calc(std::vector<ftype> NewTimeArray)
{



}

void InputTable::imped_calc(std::vector<ftype> NewFrequencyArray)
{

}

