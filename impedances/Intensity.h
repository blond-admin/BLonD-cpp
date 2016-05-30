/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INTENSITY_H_
#define IMPEDANCES_INTENSIY_H_

//class Intensity;

#include "configuration.h"
#include <complex>
#include <vector>


class Intensity {
public:
   //  *Time array of the wake in [s]*
   std::vector<ftype> fTimeArray;
   //  *Frequency array of the impedance in [Hz]*
   std::vector<ftype> fFreqArray;
   //  *Wake array in* [:math:`\Omega / s`]
   std::vector<ftype> fWake;
   //  *Impedance array in* [:math:`\Omega`]
   std::vector<complex_t> fImpedance;

   Intensity() {};
   virtual void wake_calc(std::vector<ftype> &NewTimeArray) = 0;
   virtual void imped_calc(std::vector<ftype> &NewFrequencyArray) = 0;
   virtual ~Intensity() {};
};

class Resonators: public Intensity {
public:
   // *Shunt impepdance in* [:math:`\Omega`]
   std::vector<ftype> fRS;
   // *Resonant frequency in [Hz]*
   std::vector<ftype> fFrequencyR;
   //  *Resonant angular frequency in [rad/s]*
   std::vector<ftype> fOmegaR;
   //  *Quality factor*
   std::vector<ftype> fQ;
   unsigned int fNResonators;

   void wake_calc(std::vector<ftype> &NewTimeArray);
   void imped_calc(std::vector<ftype> &NewFrequencyArray);
   Resonators(std::vector<ftype> &RS,
              std::vector<ftype> &FrequencyR, std::vector<ftype> &Q);
   ~Resonators() {} ;
};


class InputTable: public Intensity {
public:
   std::vector<ftype> fFrequencyArrayLoaded;
   std::vector<ftype> fReZArrayLoaded;
   std::vector<ftype> fImZArrayLoaded;
   std::vector<complex_t> fImpedanceLoaded;
   std::vector<ftype> fWakeArray;

   void wake_calc(std::vector<ftype> &NewTimeArray);
   void imped_calc(std::vector<ftype> &NewFrequencyArray);
   InputTable(std::vector<ftype> &input1, std::vector<ftype> &input2,
                       std::vector<ftype> input3 = std::vector<ftype>());
   //InputTable(std::vector<ftype> &input1, std::vector<ftype> &input2) {
   //   std::vector<ftype> v;
   //   InputTable(input1, input2, v);
   //};
   ~InputTable() {};
};

#endif /* IMPEDANCES_INTENSITY_H_ */