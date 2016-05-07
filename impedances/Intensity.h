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

class Intensity {
protected:
   //  *Time array of the wake in [s]*
   std::vector<ftype> fTimeArray;
   //  *Frequency array of the impedance in [Hz]*
   std::vector<ftype> fFreqArray;
   //  *Wake array in* [:math:`\Omega / s`]
   std::vector<ftype> fWake;
   //  *Impedance array in* [:math:`\Omega`]
   std::vector<ftype> fImpedance;

public:
   Intensity() {};
   virtual ~Intensity {};
   virtual void wake_calc(std::vector<ftype> NewTimeArray);
   virtual void imped_calc(std::vector<ftype> NewFrequencyArray);

};

class Resonators: public Intensity {
private:
   
   // *Shunt impepdance in* [:math:`\Omega`]
   std::vector<ftype> fRS;
   // *Resonant frequency in [Hz]*
   std::vector<ftype> fFrequencyR;
   //  *Resonant angular frequency in [rad/s]*
   std::vector<ftype> fOmegaR;
   //  *Quality factor*
   std::vector<ftype> fQ;
   uint fNResonators;

public:
   void wake_calc(std::vector<ftype> NewTimeArray);
   void imped_calc(std::vector<ftype> NewFrequencyArray);
   Resonators(int NResonators, std::vector<ftype> RS, std::vector<ftype> FrequencyR, std::vector<ftype> Q));
   ~Resonators();
};


class InputTable: public Intensity {
private:
   std::vector<ftype> fFrequencyArrayLoaded;
   std::vector<ftype> fReZArrayLoaded;
   std::vector<ftype> fImZArrayLoaded;
   std::vector<ftype> fImpedanceLoaded;

public:
   void wake_calc(std::vector<ftype> NewTimeArray);
   void imped_calc(std::vector<ftype> NewFrequencyArray);
   InputTable();
   ~InputTable();
};

#endif /* IMPEDANCES_INTENSITY_H_ */