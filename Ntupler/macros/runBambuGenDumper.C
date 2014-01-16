//--------------------------------------------------------------------------------------------------
// 
// ROOT macro to print generator level information from BAMBU
//
//==================================================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "EWKAna/Ntupler/interface/BambuGenDumperMod.hh"
#endif

using namespace mithep;

//--------------------------------------------------------------------------------------------------
void runBambuGenDumper(const char *file = "/castor/cern.ch/user/p/paus/filefi/020/p11-h160ww2l-gf-v1g1-pu/068A6ABB-F350-E011-9DCF-00A0D1EEDDA8.root")
{  
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  BambuGenDumperMod *mod = new BambuGenDumperMod;

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(mod);
  ana->AddFile(file);
  ana->SetUseHLT(kFALSE);
  //ana->SetProcessNEvents(1);
  
  // run the analysis after successful initialisation
  ana->Run(kTRUE);
}
