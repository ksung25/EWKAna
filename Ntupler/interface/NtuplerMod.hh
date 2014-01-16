#ifndef EWKANA_NTUPLER_NTUPLERMOD_H
#define EWKANA_NTUPLER_NTUPLERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/ElectronMomentumCorrection.h"

#include "EWKAnaDefs.hh"
#include "TEventInfo.hh"
#include "TGenInfo.hh"
#include <TClonesArray.h>
#include "TMuon.hh"
#include "TElectron.hh"
#include "TDielectron.hh"
#include "TJet.hh"
#include "TPhoton.hh"
#include "TVertex.hh"

#include <vector>

class TTree;
class TFile;
class TString;

namespace mithep
{  
  class NtuplerMod : public BaseMod
  {    
    public:
      NtuplerMod(const char *name="NtuplerMod", const char *title="BAMBU to ntuple");
      ~NtuplerMod();	    
      
      void SetOutputName(const char *f)         { fOutputName = f; }
      
      void SetIsData(const Bool_t flag)         { fIsData = flag; }
      void SetUseGen(const Int_t useGen)        { fUseGen = useGen; }
      void SetSkipIfHLTFail(const Bool_t flag)  { fSkipIfHLTFail = flag; } 
      void setIs2011(const bool flag=true)      { fIs2011=flag; } 
          
      void SetMuPtMin(const Double_t pt)        { fMuPtMin = pt; }
      void SetMuPtMax(const Double_t pt)        { fMuPtMax = pt; }
      void SetMuEtaMin(const Double_t eta)      { fMuEtaMin = eta; }
      void SetMuEtaMax(const Double_t eta)      { fMuEtaMax = eta; }
      void SetEleEtMin(const Double_t pt)       { fEleEtMin = pt; }
      void SetEleEtMax(const Double_t pt)       { fEleEtMax = pt; }
      void SetEleEtaMin(const Double_t eta)     { fEleEtaMin = eta; }
      void SetEleEtaMax(const Double_t eta)     { fEleEtaMax = eta; }
      void SetMassMin(const Double_t m)         { fMassMin = m; }
      void SetMassMax(const Double_t m)         { fMassMax = m; }
      void SetJetPtMin(const Double_t et)       { fJetPtMin = et; }
      void SetPhotonEtMin(const Double_t et)    { fPhotonEtMin = et; }
      void SetMinNTracksFit(const UInt_t ntrks) { fMinNTracksFit = ntrks; }
      void SetMinNdof(const Double_t ndof)      { fMinNdof = ndof; }
      void SetMaxAbsZ(const Double_t z)         { fMaxAbsZ = z; }
      void SetMaxRho(const Double_t rho)        { fMaxRho = rho; }
      void SetPrintHLT(const Bool_t flag)       { fPrintTable = flag; }
      void SetFSRMode(const Int_t mode)         { fFSRMode = mode; }
      
      void AddTrigger(const char* name, const ULong_t id,
                      const char* objName1="", const ULong_t objId1=0, const Double_t minPt1=0,
		      const char* objName2="", const ULong_t objId2=0, const Double_t minPt2=0) {
        fTriggerNamesv.push_back(name);
	fTriggerIdsv.push_back(id);
	
	fTriggerObjNames1v.push_back(objName1);
        fTriggerObjIds1v.push_back(objId1);
	fTriggerObjMinPt1v.push_back(minPt1);
	
	fTriggerObjNames2v.push_back(objName2);
        fTriggerObjIds2v.push_back(objId2);
	fTriggerObjMinPt2v.push_back(minPt2);
      }

      void AddJetCorr(const char* name) {
        fJetCorrParsv.push_back(name);
      }
            
    protected:
      void Begin();
      void BeginRun();
      void EndRun();
      void SlaveBegin();
      void SlaveTerminate();
      void Terminate();
      void Process();

      // Fill generator info data object
      void FillGenInfo(const MCParticle *p, const MCParticle *lepton1, const MCParticle *lepton2,
                       const MCParticle *pho, const UInt_t npho);
      
      // Fill muon data object
      void FillMuon(const Muon *mu);
      
      // Fill electron data object
      void FillElectron(const Electron *ele, const double scale);

      // Fill dielectron data object
      void FillDielectron(const Electron *ele1, const double scale1, const Electron *ele2, const double scale2);
            
      // Fill jet data object
      void FillJet(const PFJet *jet);
      
      // Fill photon data object
      void FillPhoton(const Photon *pho);
      void FillPhoton(const SuperCluster *sc);
      
      // Fill vertex data object
      void FillPV(const Vertex *pv);
      
      // Match muon to HLT primitive
      ULong_t MatchHLT(const Double_t eta, const Double_t phi);
      ULong_t MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi);
      
      // Check for conversion with MVF
      Bool_t IsConversion(const Electron *ele);
      
      // Compute PF isolation
      void separatePileUp(const Bool_t checkClosestZVertex);
      void computePFMuonIso(const Muon *muon, Double_t *chIso, Double_t *gammaIso, Double_t *neuHadIso, Double_t *puIso);
      void computePFEleIso(const Electron *electron, Double_t *chIso, Double_t *gammaIso, Double_t *neuHadIso);

      // Convert particle eta to detector eta (for GEN-SC matching)
      Double_t EcalEta(const Double_t eta, const Double_t z, const Double_t rho);

      // Electron momentum correction setters
      void setEleCorrParams(const TString                                    regFilename,
                            const ElectronEnergyRegression::RegressionType   regType,
                            const ElectronEnergySmearingScaling::DatasetType dataset,
                            const int                                        corrType,
                            const TString                                    scaleCorrFilename, 
                            const TString                                    smearCorrType1Filename,
                            const TString                                    smearCorrType2Filename,
                            const TString                                    smearCorrType3Filename,
                            const TString                                    linearityFilename,
                            const bool                                       doRand,
                            const int                                        seed=0xDEADBEEF,
                            const double                                     lumiRatio=1)
      {         
        fEleRegFilename            = regFilename;
        fEleRegType                = regType;
        fEleScaleSmearDataset      = dataset;
        fEleScaleSmearType         = corrType;
        fEleScaleCorrFilename      = scaleCorrFilename;
        fEleSmearCorrType1Filename = smearCorrType1Filename;
        fEleSmearCorrType2Filename = smearCorrType2Filename;
        fEleSmearCorrType3Filename = smearCorrType3Filename;
        fEleLinearityFilename      = linearityFilename;
        fEleDoRand                 = doRand;
        fEleCorrRandSeed           = seed;
        fEleCorrLumiRatio          = lumiRatio;
      }

      
      TFile                        *fOutputFile;           // output file handle
      TString                       fOutputName;           // output file name
      
      TString                       fPartName;             // MC particle collection name
      TString                       fMCEvtInfoName;        // MC event info name
      TString                       fMuonName;             // muon collection name
      TString                       fElectronName;         // electron collection name
      TString                       fPrimVtxName;          // primary vertex collection name
      TString                       fBeamSpotName;         // pointer to beam spot branch
      TString                       fPFJetName;            // particle flow jet collection name
      TString                       fPhotonName;           // photon collection name
      TString                       fTrigMaskName;         // trigger mask name
      TString                       fPFMetName;            // particle flow MET collection name
      TString                       fConversionName;       // conversion collection name
      TString                       fBarrelSCName;         // barrel superscluster collection name
      TString                       fEndcapSCName;         // endcap supercluster collection name
      TString                       fPileupName;           // pile-up info name
      TString                       fPUEnergyDensityName;  // Fastjet correction info name
      TString                       fPFCandidateName;      // particle flow candidates collection name 
      TString                       fTracksName;           // track collection name           
      
      const MCParticleCol          *fParticles;       // MC particle collection handle
      const MCEventInfo            *fMCEvtInfo;       // MC event info handle
      const MuonCol                *fMuons;           // muon collection handle
      const ElectronCol            *fElectrons;       // electron collection handle
      const VertexCol              *fPrimVerts;       // primary vertex collection handle
      const BeamSpotCol            *fBeamSpot;        // pointer to beam spot branch
      const PFJetCol               *fPFJets;          // particle flow jet collection handle
      const PhotonCol              *fPhotons;         // photon collection handle
      const TriggerMask            *fTrigMask;        // trigger mask handle
      const PFMetCol               *fPFMet;           // particle flow MET handle
      const DecayParticleCol       *fConversions;     // conversion collection handle
      const SuperClusterCol        *fBarrelSC;        // barrel supercluster collection handle
      const SuperClusterCol        *fEndcapSC;        // endcap supercluster collection handle
      const PileupInfoCol          *fPileup;          // pile-up info handle
      const PileupEnergyDensityCol *fPUEnergyDensity; // Fastjet correction info handle
      const PFCandidateCol         *fPFCandidates;    // particle flow candidates collection handle
      PFCandidateOArr              *fPFPileUp;
      PFCandidateOArr              *fPFNoPileUp;
      const TrackCol               *fTracks;          // tracks collection handle       
      
      Bool_t                  fIsData;          // flag to indicate if processing collision data
      Int_t                   fUseGen;          // flag whether to look at generator info
      Bool_t                  fIs2011;          // flag to indicate if processing 2011 sample
      Bool_t                  fPrintTable;      // flag whether to print out HLT table
      Bool_t                  fSkipIfHLTFail;   // flag whether to skip event processing if HLT does not accept
       
      Double_t                fMuPtMin;         // minimum reco muon pT
      Double_t                fMuPtMax;         // maximum reco muon pT
      Double_t                fMuEtaMin;        // minimum reco muon eta
      Double_t                fMuEtaMax;        // maximum reco muon eta
      Double_t                fEleEtMin;        // minimum electron supercluster ET
      Double_t                fEleEtMax;        // maximum electron supercluster ET
      Double_t                fEleEtaMin;       // minimum electron supercluster eta
      Double_t                fEleEtaMax;       // maximum electron supercluster eta
      Double_t                fMassMin;         // minimum dielectron mass
      Double_t                fMassMax;         // maximum dielectron mass
      Double_t                fJetPtMin;        // minimum jet ET
      Double_t                fPhotonEtMin;     // minimum photon supercluster ET
      UInt_t                  fMinNTracksFit;   // minimum number of tracks used for a good primary vertex
      Double_t                fMinNdof;         // minimum degrees of freedom for a good primary vertex
      Double_t                fMaxAbsZ;         // maximum z displacement for a good primary vertex
      Double_t                fMaxRho;          // maximum transverse displacement for a good primary vertex 
      
      TTree*                  fEventTree;       // event tree

      const Vertex*           fVertex;          // best primary vertex in the event
            
      TEventInfo              fEventInfo;       // general event information
      TGenInfo                fGenInfo;         // generator information
      TClonesArray           *fMuonArr;         // muon array
      TClonesArray           *fElectronArr;     // electron array
      TClonesArray           *fDielectronArr;   // dielectron array
      TClonesArray           *fPFJetArr;        // particle flow jet array
      TClonesArray           *fPhotonArr;       // photon array
      TClonesArray           *fPVArr;           // valid primary vertex array
      
      vector<TString>         fTriggerNamesv;   // names of triggers we're interested in 
      vector<ULong_t>         fTriggerIdsv;     // corresponding ETriggerBit value
      
      vector<TString>         fTriggerObjNames1v;
      vector<ULong_t>         fTriggerObjIds1v;
      vector<Double_t>        fTriggerObjMinPt1v;

      vector<TString>         fTriggerObjNames2v;
      vector<ULong_t>         fTriggerObjIds2v;
      vector<Double_t>        fTriggerObjMinPt2v;

      Int_t                   fFSRMode;         // flag to indicate how to parse GEN-level tree (depends on MC tool)
                                                // 0 -> PYTHIA: traverse down tree from status=3 boson
                                                // 1 -> HORACE: look at immediate daughters of status=3 boson
                                                // 2 -> PHOTOS: look at immediate daughters of status=2 boson 
      
      vector<TString>         fJetCorrParsv;
      FactorizedJetCorrector *fJetCorrector;    // CMSSW class to handle jet corrections
      
      ElectronIDMVA          *fEleMVA;  
      
      RunLumiSet                       fRunLumiSet;
      RunLumiRangeMap::RunLumiPairType fLastRunLumi;

      // Electron momentum correction parameters
      ElectronMomentumCorrection                  fEleCorr;                    // electron momentum corrector class
      TString                                     fEleRegFilename;             // electron energy regression weights file name
      ElectronEnergyRegression::RegressionType    fEleRegType;                 // electron energy regression type
      ElectronEnergySmearingScaling::DatasetType  fEleScaleSmearDataset;       // dataset type for electron scale/smear correction
      Int_t                                       fEleScaleSmearType;          // electron scale/smear correction type
      TString                                     fEleScaleCorrFilename;       // electron scale correction file name
      TString                                     fEleSmearCorrType1Filename;  // type 1 electron smear correction file name
      TString                                     fEleSmearCorrType2Filename;  // type 2 electron smear correction file name
      TString                                     fEleSmearCorrType3Filename;  // type 3 electron smear correction file name
      TString                                     fEleLinearityFilename;       // linearity correction data file name
      Bool_t                                      fEleDoRand;                  // flag to do randomized smearing
      Int_t                                       fEleCorrRandSeed;            // seed for randomized smearing
      Double_t                                    fEleCorrLumiRatio;           // fraction of total luminosity from 2012D (for MC)
      
    ClassDef(NtuplerMod,1)
  };
}

#endif
