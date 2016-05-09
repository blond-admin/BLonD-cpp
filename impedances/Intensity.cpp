/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 15:25:33
*/

#include "Intensity.h"
#include "utilities.h"
#include "constants.h"


Resonators::Resonators(std::vector<ftype> RS, std::vector<ftype> FrequencyR,
                       std::vector<ftype> Q)
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



void Resonators::wake_calc(std::vector<ftype> NewTimeArray)
{
   /*
   * Wake calculation method as a function of time.*
   */
   fTimeArray = NewTimeArray;
   fWake.resize(fTimeArray.size());
   //std::fill_n(fWake, fNResonators, 0);

   for (uint i = 0; i < fNResonators; ++i) {
      ftype alpha = fOmegaR[i] / (2 * fQ[i]);
      ftype omega_bar = sqrt(fOmegaR[i] * fOmegaR[i] - alpha * alpha);
      for (uint j = 0; j < fWake.size(); ++j) {
         ftype temp = fTimeArray[j];
         fWake[i] += (std::signbit(temp) + 1) * fRS[i] * alpha *
                     exp(-alpha * temp) * (cos(omega_bar * temp)
                                           - alpha / omega_bar * sin(omega_bar * temp));
      }
   }

}

/*
constexpr std::complex<double> operator""i(long double d)
{
    return std::complex<double>{0.0, static_cast<double>(d)};
}
*/

void Resonators::imped_calc(std::vector<ftype> NewFrequencyArray)
{
   /*
   * Impedance calculation method as a function of frequency.*
   */
   //using namespace std::complex_literals;

   fFreqArray = NewFrequencyArray;
   fImpedance.resize(fFreqArray.size());
   fImpedance[0] = complex_t(0, 0);
   for (uint i = 1; i < fImpedance.size(); i++)
      for (uint j = 0; j < fNResonators; ++j) {
         fImpedance[i] = fImpedance[i] +
                         complex_t(fRS[j], 0) /
                         complex_t(1, fQ[j] *
                                   (fFreqArray[i] / fFrequencyR[j] -
                                    fFrequencyR[j] / fFreqArray[i]));
         //fImpedance[i] += fRS[j] / (1 + 1i* fQ[j] *
         //                           (fFreqArray[i] / fFrequencyR[j] -
         //                           fFrequencyR[j]/ fFreqArray[i]));
      }


}

InputTable::InputTable()
{

}


void InputTable::wake_calc(std::vector<ftype> NewTimeArray)
{
   return;
}

void InputTable::imped_calc(std::vector<ftype> NewFrequencyArray)
{
   return;
}

