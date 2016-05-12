/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INDUCEDVOLTAGE_H_
#define IMPEDANCES_INDUCEDVOLTAGE_H_

//class InducedVoltage;

#include "configuration.h"
//#include <complex>
#include <vector>
#include "Intensity.h"

//typedef std::complex<float> complex_t;



//unsigned int next_regular(unsigned int target);

class InducedVoltage {
public:
   ftype fInducedVoltage;

   InducedVoltage() {};
   virtual void track() = 0;
   virtual void reprocess() = 0;
   virtual void induced_voltage_generation() = 0;
   virtual ~InducedVoltage() {};
};



class InducedVoltageTime: public InducedVoltage {
public:
   enum time_or_freq {
      time, freq
   };

   std::vector<Intensity *> fWakeSourceList;
   std::vector<ftype> fTimeArray;
   std::vector<ftype> fTotalWake;
   unsigned int fCut;
   unsigned int fShape;
   time_or_freq fTimeOrFreq;


   void track();
   void sum_wakes(std::vector<ftype>& v);
   void reprocess();
   void induced_voltage_generation();
   InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
                      time_or_freq TimeOrFreq = freq);
   ~InducedVoltageTime() {};
};


class InducedVoltageFreq: public InducedVoltage {
public:
   void track();
   void sum_impedances();
   void reprocess();
   void induced_voltage_generation();
   InducedVoltageFreq();
   ~InducedVoltageFreq() {};
};

#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */