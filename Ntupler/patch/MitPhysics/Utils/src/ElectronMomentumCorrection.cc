#include "MitPhysics/Utils/interface/ElectronMomentumCorrection.h"
#include <iostream>

using namespace mithep;

//--------------------------------------------------------------------------------------------------
void ElectronMomentumCorrection::initialize(
  const char   *regressionFilename,			     // regression weights file name
  const ElectronEnergyRegression::RegressionType type,       // regression type
  const ElectronEnergySmearingScaling::DatasetType dataset,  // dataset type
  const int	corrType,  			             // correction type
  const char   *scalesFilename, 			     // scale correction
  const char   *smearsType1Filename,			     // resolution correction type 1
  const char   *smearsType2Filename,			     // resolution correction type 2
  const char   *smearsType3Filename,			     // resolution correction type 3
  const char   *linearityFilename,                           // linearity correction data
  const bool	doRand,			                     // flag to toggle randomization
  const int	seed, 		                             // seed for randomization
  const double  lumiRatio)  			             // fraction of total luminosity from 2012D
{
  fRegression.initialize(regressionFilename, type);
  fSmearScale.initialize(dataset, corrType, scalesFilename, smearsType1Filename, smearsType2Filename, smearsType3Filename, doRand, seed, lumiRatio);
  fEpCombine.initialize(regressionFilename);
  fLinearity.initialize(dataset, linearityFilename);
  
  fIsInitialized = fRegression.isInitialized() && fSmearScale.isInitialized() && fEpCombine.isInitialized() && fLinearity.isInitialized();
}


//--------------------------------------------------------------------------------------------------
std::pair<double,double> ElectronMomentumCorrection::ElectronMomentumCorrection::evaluate(
  const Electron     *ele,
  const double        rho,
  const int	      nvertices,
  const unsigned int  runNum,
  const bool          printDebug)
{
  if(printDebug) {
    std::cout << "[ElectronMomentumCorrection]" << std::endl;
    std::cout << " Electron pt = " << ele->Pt() << " eta = " << ele->Eta() << " phi = " << ele->Phi() << std::endl;
  }
  
  std::pair<double,double>
  result = fRegression.evaluate(ele,
                                rho,
				nvertices,
				printDebug);  
  
  result = fSmearScale.evaluate(ele,
                                result.first,
			        result.second,
			        runNum,
				printDebug);
  
  result = fEpCombine.evaluate(ele,
                               result.first,
			       result.second,
			       printDebug);
  
  result.first *= fLinearity.corrScale(ele, result.first, printDebug);
  
  return result;
}
