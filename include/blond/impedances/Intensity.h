/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INTENSITY_H_
#define IMPEDANCES_INTENSITY_H_

// class API  Intensity;

#include <blond/configuration.h>
#include <blond/utilities.h>
#include <vector>

class API Intensity {
  public:
    //  *Time array of the wake in [s]*
    f_vector_t fTimeArray;
    //  *Frequency array of the impedance in [Hz]*
    f_vector_t fFreqArray;
    //  *Wake array in* [:math:`\Omega / s`]
    f_vector_t fWake;
    //  *Impedance array in* [:math:`\Omega`]
    complex_vector_t fImpedance;

    Intensity(){};
    virtual void wake_calc(const f_vector_t& NewTimeArray) = 0;
    virtual void imped_calc(const f_vector_t& NewFrequencyArray) = 0;
    virtual ~Intensity(){};
};

class API Resonators : public Intensity {
  public:
    // *Shunt impepdance in* [:math:`\Omega`]
    f_vector_t fRS;
    // *Resonant frequency in [Hz]*
    f_vector_t fFrequencyR;
    //  *Resonant angular frequency in [rad/s]*
    f_vector_t fOmegaR;
    //  *Quality factor*
    f_vector_t fQ;
    unsigned int fNResonators;

    void wake_calc(const f_vector_t& NewTimeArray);
    void imped_calc(const f_vector_t& NewFrequencyArray);
    Resonators(f_vector_t& RS, f_vector_t& FrequencyR, f_vector_t& Q);
    ~Resonators();
};

class API InputTable : public Intensity {
  public:
    f_vector_t fFrequencyArrayLoaded;
    f_vector_t fReZArrayLoaded;
    f_vector_t fImZArrayLoaded;
    complex_vector_t fImpedanceLoaded;
    f_vector_t fWakeArray;

    void wake_calc(const f_vector_t& NewTimeArray);
    void imped_calc(const f_vector_t& NewFrequencyArray);
    InputTable(const f_vector_t& input1, const f_vector_t& input2,
               const f_vector_t input3 = f_vector_t());
    ~InputTable();
};

#endif /* IMPEDANCES_INTENSITY_H_ */
