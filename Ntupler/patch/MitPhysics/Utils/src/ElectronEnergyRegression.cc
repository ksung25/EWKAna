#include "MitPhysics/Utils/interface/ElectronEnergyRegression.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "Cintex/Cintex.h"
#include <TFile.h>
#include <cmath>
#include <cassert>
#include <iostream>

ClassImp(mithep::ElectronEnergyRegression)

using namespace mithep;


ElectronEnergyRegression::ElectronEnergyRegression() : 
  fIsInitialized(kFALSE),
  fVersionType(kNoTrkVar),
  forestCorrection_eb(0), 
  forestCorrection_ee(0), 
  forestUncertainty_eb(0), 
  forestUncertainty_ee(0) {
}


ElectronEnergyRegression::~ElectronEnergyRegression()
{
  if(forestCorrection_eb)  delete forestCorrection_eb;
  if(forestCorrection_ee)  delete forestCorrection_ee;
  if(forestUncertainty_eb) delete forestUncertainty_eb;
  if(forestUncertainty_ee) delete forestUncertainty_ee;
}


void ElectronEnergyRegression::initialize(const char *regressionFilename, RegressionType type) {

  ROOT::Cintex::Cintex::Enable();  // Needed to read non-TObject classes from a ROOT file!
  
  // Load "forests"
  TFile file(regressionFilename);
  forestCorrection_eb  = (GBRForest*)file.Get("EBCorrection");  assert(forestCorrection_eb);
  forestCorrection_ee  = (GBRForest*)file.Get("EECorrection");  assert(forestCorrection_ee);
  forestUncertainty_eb = (GBRForest*)file.Get("EBUncertainty"); assert(forestUncertainty_eb);
  forestUncertainty_ee = (GBRForest*)file.Get("EEUncertainty"); assert(forestUncertainty_ee);  
  
  // Updating type and marking as initialized
  fVersionType   = type;
  fIsInitialized = true;
}

std::pair<double,double> ElectronEnergyRegression::evaluate(const Electron *ele, double rho, int nvertices, bool printDebug)
{  
  assert(ele);
  
  // Checking if instance has been initialized
  assert(fIsInitialized);  
  
  const SuperCluster *sc   = ele->SCluster();
  const BasicCluster *seed = ele->SCluster()->Seed();
  
  //
  // Base variables
  //
  double SCRawEnergy      = sc->RawEnergy();                        // SuperCluster raw energy
  double scEta            = sc->Eta();                              // SuperCluster eta
  double scPhi            = sc->Phi();                              // SuperCluster phi
  double R9               = sc->R9();                               // SuperCluster R9
  double etawidth         = sc->EtaWidth();                         // SuperCluster width in eta
  double phiwidth         = sc->PhiWidth();                         // SuperCluster width in phi
  double NClusters        = sc->NClusters();                        // number of basic clusters in the SuperCluster
  double HoE              = ele->HadronicOverEm();                  // electron H/E
  double EtaSeed          = seed->Eta();                            // seed cluster eta
  double PhiSeed          = seed->Phi();                            // seed cluster phi
  double ESeed            = seed->Energy();                         // seed cluster energy
  double E3x3Seed         = seed->E3x3();                           // seed cluster E3x3
  double E5x5Seed         = seed->E5x5();                           // seed cluster E5x5
  double see              = sqrt(seed->CoviEtaiEta());              // SuperCluster sigma_ieta_ieta
  double spp              = sqrt(seed->CoviPhiiPhi());              // SuperCluster sigma_iphi_iphi
  double sep              = seed->CoviEtaiPhi()/see/spp;            // SuperCluster cov(ieta,iphi)/(sigma_ieta_ieta * sigma_iphi_iphi)
  double EMaxSeed         = seed->EMax();                           // seed cluster EMax
  double E2ndSeed         = seed->E2nd();                           // seed cluster E2nd
  double ETopSeed         = seed->ETop();                           // seed cluster ETop
  double EBottomSeed      = seed->EBottom();                        // seed cluster EBottom
  double ELeftSeed        = seed->ELeft();                          // seed cluster ELeft
  double ERightSeed       = seed->ERight();                         // seed cluster ERight
  double E2x5MaxSeed      = seed->E2x5Max();                        // seed cluster E2x5Max
  double E2x5TopSeed      = seed->E2x5Top();                        // seed cluster E2x5Top
  double E2x5BottomSeed   = seed->E2x5Bottom();                     // seed cluster E2x5Bottom
  double E2x5LeftSeed     = seed->E2x5Left();                       // seed cluster E2x5Left
  double E2x5RightSeed    = seed->E2x5Right();                      // seed cluster E2x5Right
  double IEtaSeed         = seed->IEta();                           // seed cluster IEta
  double IPhiSeed         = seed->IPhi();                           // seed cluster IPhi
  double EtaCrySeed       = seed->EtaCry();                         // eta of highest energy crystal in seed cluster
  double PhiCrySeed       = seed->PhiCry();                         // phi of highest energy crystal in seed cluster
  double PreShowerOverRaw = sc->PreshowerEnergy()/sc->RawEnergy();  // (SuperCluster PreShower energy) / (SuperCluster raw energy)  
  int	 IsEcalDriven     = ele->IsEcalDriven();                    // IsEcalDriven electron?
  
  //
  // Track variables
  //
  double pt                     = ele->Pt();                                              // electron pT
  double GsfTrackPIn            = ele->PIn();                                             // electron GSF track momentum at PCA to beam line
  double fbrem                  = fmax(ele->FBrem(), -1.0);                               // electron fbrem
  double Charge                 = ele->Charge();                                          // electron charge
  double EoP                    = fmin(ele->ESuperClusterOverP(), 20.);                   // E/p  
  double TrackMomentumError     = ele->TrackMomentumError();                              // electron track momentum error
  double EcalEnergyError        = ele->EcalEnergyError();                                 // electron ECAL energy error
  int	 Classification         = ele->Classification();                                  // electron classification  
  double detaIn                 = fmin(fabs(ele->DeltaEtaSuperClusterTrackAtVtx()),0.6);  // Delta-eta between SC and track position at PCA to beam line
  double dphiIn                 = ele->DeltaPhiSuperClusterTrackAtVtx();                  // Delta-phi between SC and track position at PCA to beam line
  double detaCalo               = ele->DeltaEtaSeedClusterTrackAtCalo();                  // Delta-eta between SC and track position at calorimeter extracted from outermost state
  double dphiCalo               = ele->DeltaPhiSeedClusterTrackAtCalo();                  // Delta-phi between SC and track position at calorimeter extracted from outermost state
  double GsfTrackChiSqr         = (ele->GsfTrk()->Chi2())/(ele->GsfTrk()->Ndof());        // GSF track fit normalized chisquare (chi2/ndf)
  double ElectronEnergyOverPout = fmin(ele->EEleClusterOverPout(), 20.);                  // (electron cluster energy)/(track momentum at calorimeter)
  
  double KFTrackNLayers = 0;  // number of layers with measurements in CTF track
  if(ele->HasTrackerTrk()) {
    const Track *trk = ele->TrackerTrk();
    
    // loop through layers and count
    if(trk->Hit(Track::PXB1)) KFTrackNLayers++; 
    if(trk->Hit(Track::PXB2)) KFTrackNLayers++; 
    if(trk->Hit(Track::PXB3)) KFTrackNLayers++; 
    
    if(trk->Hit(Track::PXF1)) KFTrackNLayers++; 
    if(trk->Hit(Track::PXF2)) KFTrackNLayers++; 
    
    if(trk->Hit(Track::TIB1) || trk->Hit(Track::TIB1S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TIB2) || trk->Hit(Track::TIB2S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TIB3))                           KFTrackNLayers++; 
    if(trk->Hit(Track::TIB4))                           KFTrackNLayers++; 
    
    if(trk->Hit(Track::TID1) || trk->Hit(Track::TID1S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TID2) || trk->Hit(Track::TID2S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TID3) || trk->Hit(Track::TID3S)) KFTrackNLayers++; 
    
    if(trk->Hit(Track::TOB1) || trk->Hit(Track::TOB1S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TOB2) || trk->Hit(Track::TOB2S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TOB3))                           KFTrackNLayers++; 
    if(trk->Hit(Track::TOB4))                           KFTrackNLayers++; 
    if(trk->Hit(Track::TOB5))                           KFTrackNLayers++; 
    if(trk->Hit(Track::TOB6))                           KFTrackNLayers++; 
    
    if(trk->Hit(Track::TEC1) || trk->Hit(Track::TEC1S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC2) || trk->Hit(Track::TEC2S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC3) || trk->Hit(Track::TEC3S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC4) || trk->Hit(Track::TEC4S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC5) || trk->Hit(Track::TEC5S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC6) || trk->Hit(Track::TEC6S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC7) || trk->Hit(Track::TEC7S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC8) || trk->Hit(Track::TEC8S)) KFTrackNLayers++; 
    if(trk->Hit(Track::TEC9) || trk->Hit(Track::TEC9S)) KFTrackNLayers++; 
 
  } else {
    KFTrackNLayers = -1;
  }
    
  //
  // Sub-cluster variables
  //
  double ESubs=0;  // Energy sum of all basic clusters except seed
  for(UInt_t iclus=0; iclus<sc->NClusters(); iclus++) {
    if(sc->Cluster(iclus) == seed) continue;
    ESubs+=sc->Cluster(iclus)->Energy();
  }    

  // sub-clusters after seed cluster (assumes clusters are sorted by energy)
  const BasicCluster *sub1=0, *sub2=0, *sub3=0;
  for(UInt_t iclus=0; iclus<sc->NClusters(); iclus++) {
    if(sc->Cluster(iclus) == seed) continue;
    if     (sub1==0) sub1 = sc->Cluster(iclus);
    else if(sub2==0) sub2 = sc->Cluster(iclus);
    else if(sub3==0) sub3 = sc->Cluster(iclus);
  }

  double EPshwSubs=0;       // Energy sum of preshower clusters
  for(UInt_t iclus=0; iclus<sc->NPsClusts(); iclus++) {
    EPshwSubs+=sc->PsClust(iclus)->Energy();
  }
  
  // preshower clusters (assumes clusters are sorted by energy)
  const PsCluster *pshwsub1 = (sc->NPsClusts()>0) ? sc->PsClust(0) : 0;
  const PsCluster *pshwsub2 = (sc->NPsClusts()>1) ? sc->PsClust(1) : 0;
  const PsCluster *pshwsub3 = (sc->NPsClusts()>2) ? sc->PsClust(2) : 0;
  
  double isEtaGap      = ele->IsEBEtaGap();                     // is in EB eta module gap
  double isPhiGap      = ele->IsEBPhiGap();                     // is in EB phi module gap
  double isDeeGap      = ele->IsEEDeeGap();                     // is in EE Dee gap
  double ESub1         = sub1 ? sub1->Energy() : 0.;            // 1st sub-cluster energy
  double EtaSub1       = sub1 ? sub1->Eta()    : 999.;          // 1st sub-cluster eta
  double PhiSub1       = sub1 ? sub1->Phi()    : 999.;          // 1st sub-cluster phi
  double EMaxSub1      = sub1 ? sub1->EMax()   : 0.;            // 1st sub-cluster EMax
  double E3x3Sub1      = sub1 ? sub1->E3x3()   : 0.;            // 1st sub-cluster E3x3
  double ESub2         = sub2 ? sub2->Energy() : 0.;            // 2nd sub-cluster energy
  double EtaSub2       = sub2 ? sub2->Eta()    : 999.;          // 2nd sub-cluster eta
  double PhiSub2       = sub2 ? sub2->Phi()    : 999.;          // 2nd sub-cluster phi
  double EMaxSub2      = sub2 ? sub2->EMax()   : 0.;            // 2nd sub-cluster EMax
  double E3x3Sub2      = sub2 ? sub2->E3x3()   : 0.;            // 2nd sub-cluster E3x3
  double ESub3         = sub3 ? sub3->Energy() : 0.;            // 3rd sub-cluster energy
  double EtaSub3       = sub3 ? sub3->Eta()    : 999.;          // 3rd sub-cluster eta
  double PhiSub3       = sub3 ? sub3->Phi()    : 999.;          // 3rd sub-cluster phi
  double EMaxSub3      = sub3 ? sub3->EMax()   : 0.;            // 3rd sub-cluster EMax
  double E3x3Sub3      = sub3 ? sub3->E3x3()   : 0.;            // 3rd sub-cluster E3x3
  double NPshwClusters = sc->NPsClusts();                       // number of preshower clusters
  double EPshwSub1     = pshwsub1 ? pshwsub1->Energy() : 0.;    // 1st preshower cluster energy
  double EtaPshwSub1   = pshwsub1 ? pshwsub1->Eta()    : 999.;  // 1st preshower cluster eta
  double PhiPshwSub1   = pshwsub1 ? pshwsub1->Phi()    : 999.;  // 1st preshower cluster phi
  double EPshwSub2     = pshwsub2 ? pshwsub2->Energy() : 0.;    // 2nd preshower cluster energy
  double EtaPshwSub2   = pshwsub2 ? pshwsub2->Eta()    : 999.;  // 2nd preshower cluster eta
  double PhiPshwSub2   = pshwsub2 ? pshwsub2->Phi()    : 999.;  // 2nd preshower cluster phi
  double EPshwSub3     = pshwsub3 ? pshwsub3->Energy() : 0.;    // 3rd preshower cluster energy
  double EtaPshwSub3   = pshwsub3 ? pshwsub3->Eta()    : 999.;  // 3rd preshower cluster eta
  double PhiPshwSub3   = pshwsub3 ? pshwsub3->Phi()    : 999.;  // 3rd preshower cluster phi
  
    
  bool isBarrel = false;
  uint nvars = 0;
  float *vals = 0; 
  
  if(fVersionType == kNoTrkVar) {  ///// Regression without tracker variables /////    
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 38 : 31;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IEtaSeed;
      vals[31] = IPhiSeed;
      vals[32] = ((int) IEtaSeed)%5;
      vals[33] = ((int) IPhiSeed)%2;
      vals[34] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[35] = ((int) IPhiSeed)%20;
      vals[36] = EtaCrySeed;
      vals[37] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = PreShowerOverRaw;
    }
    
  } else if(fVersionType == kNoTrkVarV1) {  ///// Regression without tracker variables (V1) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 39 : 32;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = IEtaSeed;
      vals[32] = IPhiSeed;
      vals[33] = ((int) IEtaSeed)%5;
      vals[34] = ((int) IPhiSeed)%2;
      vals[35] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[36] = ((int) IPhiSeed)%20;
      vals[37] = EtaCrySeed;
      vals[38] = PhiCrySeed;    
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = PreShowerOverRaw;    
    }
  
  } else if(fVersionType == kWithTrkVar) {  ///// Regression with tracker variables /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 43 : 36;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = pt;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = IEtaSeed;
      vals[36] = IPhiSeed;
      vals[37] = ((int) IEtaSeed)%5;
      vals[38] = ((int) IPhiSeed)%2;
      vals[39] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[40] = ((int) IPhiSeed)%20;
      vals[41] = EtaCrySeed;
      vals[42] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = pt;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = PreShowerOverRaw;    
    }
  
  } else if(fVersionType == kWithTrkVarV1) {  ///// Regression with tracker variables (V1) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 46 : 39;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = IEtaSeed;
      vals[39] = IPhiSeed;
      vals[40] = ((int) IEtaSeed)%5;
      vals[41] = ((int) IPhiSeed)%2;
      vals[42] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[43] = ((int) IPhiSeed)%20;
      vals[44] = EtaCrySeed;
      vals[45] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = PreShowerOverRaw;
    }
  
  } else if(fVersionType == kWithTrkVarV2) {  ///// Regression with tracker variables (V2) /////
    isBarrel = (fabs(scEta) <= 1.479);
    nvars = isBarrel ? 53 : 46;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = detaIn;
      vals[39] = dphiIn;
      vals[40] = detaCalo;
      vals[41] = dphiCalo;
      vals[42] = GsfTrackChiSqr;
      vals[43] = KFTrackNLayers;
      vals[44] = ElectronEnergyOverPout;
      vals[45] = IEtaSeed;
      vals[46] = IPhiSeed;
      vals[47] = ((int) IEtaSeed)%5;
      vals[48] = ((int) IPhiSeed)%2;
      vals[49] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[50] = ((int) IPhiSeed)%20;
      vals[51] = EtaCrySeed;
      vals[52] = PhiCrySeed;
    
    } else {
      vals[0]  = SCRawEnergy;
      vals[1]  = scEta;
      vals[2]  = scPhi;
      vals[3]  = R9;
      vals[4]  = E5x5Seed/SCRawEnergy;
      vals[5]  = etawidth;
      vals[6]  = phiwidth;
      vals[7]  = NClusters;
      vals[8]  = HoE;
      vals[9]  = rho;
      vals[10] = nvertices;
      vals[11] = EtaSeed - scEta;
      vals[12] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[13] = ESeed/SCRawEnergy;
      vals[14] = E3x3Seed/ESeed;
      vals[15] = E5x5Seed/ESeed;
      vals[16] = see;
      vals[17] = spp;
      vals[18] = sep;
      vals[19] = EMaxSeed/ESeed;
      vals[20] = E2ndSeed/ESeed;
      vals[21] = ETopSeed/ESeed;
      vals[22] = EBottomSeed/ESeed;
      vals[23] = ELeftSeed/ESeed;
      vals[24] = ERightSeed/ESeed;
      vals[25] = E2x5MaxSeed/ESeed;
      vals[26] = E2x5TopSeed/ESeed;
      vals[27] = E2x5BottomSeed/ESeed;
      vals[28] = E2x5LeftSeed/ESeed;
      vals[29] = E2x5RightSeed/ESeed;
      vals[30] = IsEcalDriven;
      vals[31] = GsfTrackPIn;
      vals[32] = fbrem;
      vals[33] = Charge;
      vals[34] = EoP;
      vals[35] = TrackMomentumError/GsfTrackPIn;
      vals[36] = EcalEnergyError/SCRawEnergy;
      vals[37] = Classification;
      vals[38] = detaIn;
      vals[39] = dphiIn;
      vals[40] = detaCalo;
      vals[41] = dphiCalo;
      vals[42] = GsfTrackChiSqr;
      vals[43] = KFTrackNLayers;
      vals[44] = ElectronEnergyOverPout;
      vals[45] = PreShowerOverRaw;    
    } 
  
  } else if(fVersionType == kWithSubCluVar) {  ///// Regression with subcluster variables and without track variables /////
    isBarrel = ele->IsEB();
    nvars = isBarrel ? 61 : 65;
    vals = new float[nvars];
    
    if(isBarrel) {
      vals[0]  = rho;
      vals[1]  = nvertices;
      vals[2]  = IsEcalDriven;
      vals[3]  = isEtaGap;
      vals[4]  = isPhiGap;
      vals[5]  = isDeeGap;
      vals[6]  = SCRawEnergy;
      vals[7]  = scEta;
      vals[8]  = scPhi;
      vals[9]  = R9;
      vals[10] = etawidth;
      vals[11] = phiwidth;
      vals[12] = NClusters;
      vals[13] = HoE;
      vals[14] = EtaSeed - scEta;
      vals[15] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[16] = ESeed/SCRawEnergy;
      vals[17] = E3x3Seed/ESeed;
      vals[18] = E5x5Seed/SCRawEnergy;
      vals[19] = E5x5Seed/ESeed;
      vals[20] = EMaxSeed/ESeed;
      vals[21] = E2ndSeed/ESeed;
      vals[22] = ETopSeed/ESeed;
      vals[23] = EBottomSeed/ESeed;
      vals[24] = ELeftSeed/ESeed;
      vals[25] = ERightSeed/ESeed;
      vals[26] = E2x5MaxSeed/ESeed;
      vals[27] = E2x5TopSeed/ESeed;
      vals[28] = E2x5BottomSeed/ESeed;
      vals[29] = E2x5LeftSeed/ESeed;
      vals[30] = E2x5RightSeed/ESeed;
      vals[31] = see;
      vals[32] = spp;
      vals[33] = sep;
      vals[34] = phiwidth/etawidth;
      vals[35] = (ELeftSeed+ERightSeed==0. ? 0. : (ELeftSeed-ERightSeed)/(ELeftSeed+ERightSeed));
      vals[36] = (ETopSeed+EBottomSeed==0. ? 0. : (ETopSeed-EBottomSeed)/(ETopSeed+EBottomSeed));
      vals[37] = ESubs/SCRawEnergy;
      vals[38] = ESub1/SCRawEnergy;
      vals[39] = (NClusters<=1 ? 999. : EtaSub1-EtaSeed);
      vals[40] = (NClusters<=1 ? 999. : atan2(sin(PhiSub1-PhiSeed),cos(PhiSub1-PhiSeed)));
      vals[41] = (NClusters<=1 ? 0.   : EMaxSub1/ESub1);
      vals[42] = (NClusters<=1 ? 0.   : E3x3Sub1/ESub1);
      vals[43] = ESub2/SCRawEnergy;
      vals[44] = (NClusters<=2 ? 999. : EtaSub2-EtaSeed);
      vals[45] = (NClusters<=2 ? 999. : atan2(sin(PhiSub2-PhiSeed),cos(PhiSub2-PhiSeed)));
      vals[46] = (NClusters<=2 ? 0.   : EMaxSub2/ESub2);
      vals[47] = (NClusters<=2 ? 0.   : E3x3Sub2/ESub2);
      vals[48] = ESub3/SCRawEnergy;
      vals[49] = (NClusters<=3 ? 999. : EtaSub3-EtaSeed);
      vals[50] = (NClusters<=3 ? 999. : atan2(sin(PhiSub3-PhiSeed),cos(PhiSub3-PhiSeed)));
      vals[51] = (NClusters<=3 ? 0.   : EMaxSub3/ESub3);
      vals[52] = (NClusters<=3 ? 0.   : E3x3Sub3/ESub3);
      vals[53] = IEtaSeed;
      vals[54] = ((int) IEtaSeed)%5;
      vals[55] = (abs(IEtaSeed)<=25)*(((int)IEtaSeed)%25) + (abs(IEtaSeed)>25)*(((int) (IEtaSeed-25*abs(IEtaSeed)/IEtaSeed))%20);
      vals[56] = IPhiSeed;
      vals[57] = ((int) IPhiSeed)%2;
      vals[58] = ((int) IPhiSeed)%20;
      vals[59] = EtaCrySeed;
      vals[60] = PhiCrySeed;
    
    } else {
      vals[0]  = rho;
      vals[1]  = nvertices;
      vals[2]  = IsEcalDriven;
      vals[3]  = isEtaGap;
      vals[4]  = isPhiGap;
      vals[5]  = isDeeGap;
      vals[6]  = SCRawEnergy;
      vals[7]  = scEta;
      vals[8]  = scPhi;
      vals[9]  = R9;
      vals[10] = etawidth;
      vals[11] = phiwidth;
      vals[12] = NClusters;
      vals[13] = HoE;
      vals[14] = EtaSeed - scEta;
      vals[15] = atan2(sin(PhiSeed-scPhi),cos(PhiSeed-scPhi));
      vals[16] = ESeed/SCRawEnergy;
      vals[17] = E3x3Seed/ESeed;
      vals[18] = E5x5Seed/SCRawEnergy;
      vals[19] = E5x5Seed/ESeed;
      vals[20] = EMaxSeed/ESeed;
      vals[21] = E2ndSeed/ESeed;
      vals[22] = ETopSeed/ESeed;
      vals[23] = EBottomSeed/ESeed;
      vals[24] = ELeftSeed/ESeed;
      vals[25] = ERightSeed/ESeed;
      vals[26] = E2x5MaxSeed/ESeed;
      vals[27] = E2x5TopSeed/ESeed;
      vals[28] = E2x5BottomSeed/ESeed;
      vals[29] = E2x5LeftSeed/ESeed;
      vals[30] = E2x5RightSeed/ESeed;
      vals[31] = see;
      vals[32] = spp;
      vals[33] = sep;
      vals[34] = phiwidth/etawidth;
      vals[35] = (ELeftSeed+ERightSeed==0. ? 0. : (ELeftSeed-ERightSeed)/(ELeftSeed+ERightSeed));
      vals[36] = (ETopSeed+EBottomSeed==0. ? 0. : (ETopSeed-EBottomSeed)/(ETopSeed+EBottomSeed));
      vals[37] = ESubs/SCRawEnergy;
      vals[38] = ESub1/SCRawEnergy;
      vals[39] = (NClusters<=1 ? 999. : EtaSub1-EtaSeed);
      vals[40] = (NClusters<=1 ? 999. : atan2(sin(PhiSub1-PhiSeed),cos(PhiSub1-PhiSeed)));
      vals[41] = (NClusters<=1 ? 0.   : EMaxSub1/ESub1);
      vals[42] = (NClusters<=1 ? 0.   : E3x3Sub1/ESub1);
      vals[43] = ESub2/SCRawEnergy;
      vals[44] = (NClusters<=2 ? 999. : EtaSub2-EtaSeed);
      vals[45] = (NClusters<=2 ? 999. : atan2(sin(PhiSub2-PhiSeed),cos(PhiSub2-PhiSeed)));
      vals[46] = (NClusters<=2 ? 0.   : EMaxSub2/ESub2);
      vals[47] = (NClusters<=2 ? 0.   : E3x3Sub2/ESub2);
      vals[48] = ESub3/SCRawEnergy;
      vals[49] = (NClusters<=3 ? 999. : EtaSub3-EtaSeed);
      vals[50] = (NClusters<=3 ? 999. : atan2(sin(PhiSub3-PhiSeed),cos(PhiSub3-PhiSeed)));
      vals[51] = (NClusters<=3 ? 0.   : EMaxSub3/ESub3);
      vals[52] = (NClusters<=3 ? 0.   : E3x3Sub3/ESub3);
      vals[53] = PreShowerOverRaw;
      vals[54] = NPshwClusters;
      vals[55] = EPshwSubs/SCRawEnergy;
      vals[56] = EPshwSub1/SCRawEnergy;
      vals[57] = (NPshwClusters==0 ? 999. : EtaPshwSub1-EtaSeed);
      vals[58] = (NPshwClusters==0 ? 999. : atan2(sin(PhiPshwSub1-PhiSeed),cos(PhiPshwSub1-PhiSeed)));
      vals[59] = EPshwSub2/SCRawEnergy;
      vals[60] = (NPshwClusters<=1 ? 999. : EtaPshwSub2-EtaSeed);
      vals[61] = (NPshwClusters<=1 ? 999. : atan2(sin(PhiPshwSub2-PhiSeed),cos(PhiPshwSub2-PhiSeed)));
      vals[62] = EPshwSub3/SCRawEnergy;
      vals[63] = (NPshwClusters<=2 ? 999. : EtaPshwSub3-EtaSeed);
      vals[64] = (NPshwClusters<=2 ? 999. : atan2(sin(PhiPshwSub3-PhiSeed),cos(PhiPshwSub3-PhiSeed)));    
    }
  
  } else {
    std::cout << "Error: Regression VersionType " << fVersionType << " is not supported!" << std::endl;
    assert(0);    
  }

  // evaluate regression
  double value       = isBarrel ? (SCRawEnergy * forestCorrection_eb->GetResponse(vals))
                                : (SCRawEnergy * (1+PreShowerOverRaw) * forestCorrection_ee->GetResponse(vals));
  double uncertainty = isBarrel ? (SCRawEnergy * forestUncertainty_eb->GetResponse(vals))
                                : (SCRawEnergy * (1+PreShowerOverRaw) * forestUncertainty_ee->GetResponse(vals));
  Int_t BinIndex     = isBarrel ? 0 : 1;

  // print debug
  if(printDebug) {    
    std::cout << "[ElectronEnergyRegression]" << std::endl;
    std::cout << (isBarrel ? "Barrel :" : "Endcap :") << std::endl;
    for(uint v=0; v < nvars; ++v) std::cout << v << "=" << vals[v] << ", ";
    std::cout << std::endl;

    std::cout << "BinIndex : " << BinIndex << "\n";
    std::cout << "SCRawEnergy = " << SCRawEnergy << " : PreShowerOverRaw = " << PreShowerOverRaw << std::endl;
    std::cout << "regression energy = " << value << std::endl;
    std::cout << "regression energy uncertainty = " << uncertainty << std::endl;
  }
  
  // Cleaning up and returning
  delete[] vals;
  return std::pair<double,double>(value,uncertainty);
}
