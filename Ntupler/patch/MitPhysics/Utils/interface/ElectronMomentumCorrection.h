//--------------------------------------------------------------------------------------------------
//
// ElectronMomentumCorrection
//
// correct the four-momentum of the electron
//
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ELECTRONMOMENTUMCORRECTION_H
#define MITPHYSICS_UTILS_ELECTRONMOMENTUMCORRECTION_H

#include "MitAna/DataTree/interface/Electron.h"
#include "MitPhysics/Utils/interface/ElectronEnergyRegression.h"
#include "MitPhysics/Utils/interface/ElectronEnergySmearingScaling.h"
#include "MitPhysics/Utils/interface/ElectronEpCombination.h"
#include "MitPhysics/Utils/interface/ElectronLinearityCorrection.h"
#include <utility>

class TLorentzVector;

namespace mithep {

  class ElectronMomentumCorrection
  {
    public:
      ElectronMomentumCorrection():fIsInitialized(false){}
      ~ElectronMomentumCorrection(){}
      
      void initialize(
        const char   *regressionFilename,                          // regression weights file name
	const ElectronEnergyRegression::RegressionType type,       // regression type
	const ElectronEnergySmearingScaling::DatasetType dataset,  // dataset type
	const int     corrType,                                    // scale/smear correction type
	const char   *scalesFilename,                              // scale correction
        const char   *smearsType1Filename,	                   // resolution correction type 1
	const char   *smearsType2Filename,                         // resolution correction type 2
	const char   *smearsType3Filename,                         // resolution correction type 3
	const char   *linearityFilename,                           // linearity correction data
	const bool    doRand    = true,	                           // flag to toggle randomization
	const int     seed      = 0xDEADBEEF,                      // seed for randomization
	const double  lumiRatio = 1);                              // fraction of total luminosity from 2012D
      
      bool isInitialized() const {return fIsInitialized;}
      
      std::pair<double,double> evaluate(
        const Electron     *ele,                                    // electron object
	const double        rho,                                    // event energy density
	const int           nvertices,                              // number of primary vertices
	const unsigned int  runNum,                                 // run number
	const bool          printDebug=false);
          
    protected:
      bool fIsInitialized;
      
      ElectronEnergyRegression      fRegression;
      ElectronEnergySmearingScaling fSmearScale;
      ElectronEpCombination         fEpCombine;
      ElectronLinearityCorrection   fLinearity;        
  };
}
#endif
