/*
 * LHCNoiseFB.h
 *
 *  Created on: June 23, 2016
 *      Author: kiliakis
 *
 *  Feedback on phase noise amplitude for LHC controlled longitudinal emittance
 *  blow-up using noise injection through cavity controller or phase loop.
 *  The feedback compares the FWHM bunch length of the bunch to a target value
 *  and scales the phase noise to keep the targeted value.
 *  Activate the feedback either by passing it in RFSectionParameters or in
 *  the PhaseLoop object.
 *  Update the noise amplitude scaling using track().
 *  Pass the bunch pattern (occupied bucket numbers from 0...h-1) in buckets
 *  for multi-bunch simulations; the feedback uses the average bunch length.
 *
 */

#ifndef LLRF_LHCNOISEFB_H_
#define LLRF_LHCNOISEFB_H_


#include <configuration.h>
#include <functional>

class LHCNoiseFB {
private:
   constexpr static ftype cfwhm = std::sqrt(2.0 / std::log(2.0));
public:

   // Phase noise scaling factor. Initially 0
   ftype fX;

   // Target bunch length [s], 4-sigma value.
   ftype fBlTarg;

   // Measured bunch length [s], FWHM.
   ftype fBlMeas;

   // Feedback recursion scaling factor.*
   ftype fA;

   // Update feedback every n_update turns.*
   uint fNUpdate;

   // Switch to use constant or variable gain*
   bool fVariableGain;

   // Feedback gain [1/s]
   f_vector_t fG;

   // Bunch pattern for multi-bunch simulations
   f_vector_t fBunchPattern;

   f_vector_t fBlMeasBBB;

   std::function<void()> fFwhm;

   LHCNoiseFB(ftype bl_target, ftype gain = 0.1e9,
              ftype factor = 0.93, ftype update_frequency = 22500,
              bool variable_gain = true, f_vector_t bunch_pattern = f_vector_t());
   ~LHCNoiseFB() {};
   void track();
   ftype fwhm_interpolation(uint_vector_t index, ftype half_height);
   void fwhm_single_bunch();
   void fwhm_multi_bunch();
};


#endif /* LLRF_LHCNOISEFB_H_ */
