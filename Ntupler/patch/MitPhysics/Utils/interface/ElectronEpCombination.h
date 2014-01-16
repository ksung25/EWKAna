//--------------------------------------------------------------------------------------------------
//
// ElectronEnergySmearingScaling
//
// Helper Class for applying electron energy scale and resolution corrections
//
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ELECTRONEPCOMBINATION_H
#define MITPHYSICS_UTILS_ELECTRONEPCOMBINATION_H

#include "MitAna/DataTree/interface/Electron.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include <utility>

namespace mithep {

  class ElectronEpCombination {
    public:
      ElectronEpCombination():fForest(0){}
      ~ElectronEpCombination() { if(fForest) delete fForest; }
      
      void initialize(const char *regressionFilename);
      
      bool isInitialized() const {return fIsInitialized;}
      
      std::pair<double,double> evaluate(const Electron *ele,           // electron object
                                        const double    energy,        // electron energy
		                        const double    energyError,   // electron energy uncertainty
					const bool      printDebug=false);
    
    protected:
      bool fIsInitialized;
      GBRForest *fForest;
  };
}
#endif
