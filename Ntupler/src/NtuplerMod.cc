#include "EWKAna/Ntupler/interface/NtuplerMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/DecayParticle.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <vector>

using namespace mithep;

ClassImp(mithep::NtuplerMod)

NtuplerMod::NtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
  fOutputName    ("ntuple.root"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fMuonName      (Names::gkMuonBrn),
  fElectronName  (Names::gkElectronBrn),
  fPrimVtxName   (Names::gkPVBrn),
  fBeamSpotName  (Names::gkBeamSpotBrn),
  fPFJetName     (Names::gkPFJetBrn),
  fPhotonName    (Names::gkPhotonBrn),
  fTrigMaskName  (Names::gkHltBitBrn),
  fPFMetName     ("PFMet"),
  fConversionName(Names::gkMvfConversionBrn),
  fBarrelSCName  (Names::gkBarrelSuperClusterBrn),
  fEndcapSCName  (Names::gkEndcapSuperClusterBrn),
  fPileupName    (Names::gkPileupInfoBrn),
  fPUEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fTracksName    (Names::gkTrackBrn),
  fParticles     (0),
  fMCEvtInfo     (0),
  fMuons         (0),
  fElectrons     (0),
  fPrimVerts     (0),
  fBeamSpot      (0),
  fPFJets        (0),
  fPhotons       (0),  
  fTrigMask      (0),
  fPFMet         (0),
  fConversions   (0),
  fPileup        (0),
  fPUEnergyDensity(0),
  fPFCandidates  (0),
  fPFPileUp      (0),
  fPFNoPileUp    (0),
  fTracks        (0),
  fIsData        (kFALSE),
  fUseGen        (0),
  fIs2011        (true),
  fPrintTable    (kFALSE),
  fSkipIfHLTFail (kFALSE),
  fMuPtMin       (15),
  fMuPtMax       (1000),
  fMuEtaMin      (-3),
  fMuEtaMax      (3),
  fEleEtMin      (15),
  fEleEtMax      (1000),
  fEleEtaMin     (-3),
  fEleEtaMax     (3),
  fMassMin       (10),
  fMassMax       (7000),
  fJetPtMin      (15),
  fPhotonEtMin   (10),
  fMinNTracksFit (0),
  fMinNdof       (4),
  fMaxAbsZ       (24),
  fMaxRho        (2),
  fEventTree     (0),
  fFSRMode       (0),
  fJetCorrector  (0),
  fEleMVA        (0),
  fEleRegFilename           (""),
  fEleRegType               (ElectronEnergyRegression::kWithSubCluVar),   
  fEleScaleSmearDataset     (ElectronEnergySmearingScaling::kSummer12_LegacyPaper),
  fEleScaleSmearType        (1),
  fEleScaleCorrFilename     (""),
  fEleSmearCorrType1Filename(""),
  fEleSmearCorrType2Filename(""),
  fEleSmearCorrType3Filename(""),
  fEleDoRand                (true),
  fEleCorrRandSeed          (0xDEADBEEF),
  fEleCorrLumiRatio         (1)
{
  // Constructor
  
  // Don't write TObject part of the objects
  TEventInfo::Class()->IgnoreTObjectStreamer();
  TGenInfo::Class()->IgnoreTObjectStreamer();
  TMuon::Class()->IgnoreTObjectStreamer();
  TElectron::Class()->IgnoreTObjectStreamer();
  TDielectron::Class()->IgnoreTObjectStreamer();
  TJet::Class()->IgnoreTObjectStreamer();
  TPhoton::Class()->IgnoreTObjectStreamer();
  TVertex::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
NtuplerMod::~NtuplerMod()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void NtuplerMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,            fParticles); 
  ReqBranch(fMCEvtInfoName,       fMCEvtInfo);
  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fBeamSpotName,        fBeamSpot);
  ReqBranch(fPFJetName,           fPFJets);
  ReqBranch(fTrigMaskName,        fTrigMask);
  ReqBranch(fPFMetName,           fPFMet);
  ReqBranch(fPhotonName,          fPhotons);
  ReqBranch(fConversionName,      fConversions);
  ReqBranch(fBarrelSCName,        fBarrelSC);
  ReqBranch(fEndcapSCName,        fEndcapSC);
  ReqBranch(fPileupName,          fPileup);
  ReqBranch(fPUEnergyDensityName, fPUEnergyDensity);
  ReqBranch(fPFCandidateName,     fPFCandidates);
  ReqBranch(fTracksName,          fTracks);   

  // Pileup and NoPileup collections of PFCandidates
  fPFPileUp   = new PFCandidateOArr;
  fPFNoPileUp = new PFCandidateOArr;
  
  //
  // Set up arrays
  //
  fMuonArr       = new TClonesArray("mithep::TMuon");	    assert(fMuonArr);  
  fElectronArr   = new TClonesArray("mithep::TElectron");   assert(fElectronArr);
  fDielectronArr = new TClonesArray("mithep::TDielectron"); assert(fDielectronArr);
  fPFJetArr      = new TClonesArray("mithep::TJet");	    assert(fPFJetArr);
  fPhotonArr     = new TClonesArray("mithep::TPhoton");     assert(fPhotonArr);
  fPVArr         = new TClonesArray("mithep::TVertex");	    assert(fPVArr);
  
  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(!fIsData && fUseGen)
    fEventTree->Branch("Gen",&fGenInfo);
  
  fEventTree->Branch("Muon",      &fMuonArr);
  fEventTree->Branch("Electron",  &fElectronArr);
  fEventTree->Branch("Dielectron",&fDielectronArr);
  fEventTree->Branch("PFJet",     &fPFJetArr);
  fEventTree->Branch("Photon",    &fPhotonArr);
  fEventTree->Branch("PV",        &fPVArr);
  
  //
  // Set up jet corrections for PF jets
  //
  std::vector<JetCorrectorParameters> correctionParameters;
  for(UInt_t icorr=0; icorr<fJetCorrParsv.size(); icorr++)
    correctionParameters.push_back(JetCorrectorParameters(fJetCorrParsv[icorr].Data()));
    
  //initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);
  
  // initialize electron MVA
  fEleMVA = new ElectronIDMVA();
  fEleMVA->Initialize("BDTG method",
                      "Subdet0LowPt_WithIPInfo_BDTG.weights.xml",   // 0     < |eta| < 1,     10 < pT < 20
		      "Subdet1LowPt_WithIPInfo_BDTG.weights.xml",   // 1     < |eta| < 1.479, 10 < pT < 20
		      "Subdet2LowPt_WithIPInfo_BDTG.weights.xml",   // 1.479 < |eta| < 2.5,   10 < pT < 20
		      "Subdet0HighPt_WithIPInfo_BDTG.weights.xml",  // 0     < |eta| < 1,     pT > 20
		      "Subdet1HighPt_WithIPInfo_BDTG.weights.xml",  // 1     < |eta| < 1.479, pT > 20
		      "Subdet2HighPt_WithIPInfo_BDTG.weights.xml",  // 1.479 < |eta| < 2.5,   pT > 20
		      ElectronIDMVA::kWithIPInfo);

  fLastRunLumi = RunLumiRangeMap::RunLumiPairType(0,0);

  // initialize electron momentum corrections
  fEleCorr.initialize(fEleRegFilename.Data(),
                      fEleRegType,
                      fEleScaleSmearDataset,
                      fEleScaleSmearType,
                      fEleScaleCorrFilename.Data(),
                      fEleSmearCorrType1Filename.Data(),
                      fEleSmearCorrType2Filename.Data(),
                      fEleSmearCorrType3Filename.Data(),
                      fEleLinearityFilename.Data(),
                      fEleDoRand,
                      fEleCorrRandSeed,
                      fEleCorrLumiRatio);
  assert(fEleCorr.isInitialized());
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();
  
  delete fPFPileUp;
  delete fPFNoPileUp;

  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fMuonArr;
  delete fElectronArr;
  delete fDielectronArr;
  delete fPFJetArr;
  delete fPhotonArr;
  delete fPVArr;
  
  delete fJetCorrector;
  delete fEleMVA;

  //
  // Dump JSON file
  //
  TString jsonfname = fOutputName+TString(".json");
  if(fIsData)
    fRunLumiSet.DumpJSONFile(jsonfname.Data());
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::BeginRun()
{
  if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::Process()
{
  RunLumiRangeMap::RunLumiPairType rl(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec());
  if(rl!=fLastRunLumi) {
    fLastRunLumi = rl;
    fRunLumiSet.Add(rl);
  }
  
  //
  // Load branches
  //
  if(!fIsData && fUseGen) LoadBranch(fPartName);
  if(!fIsData && fUseGen) LoadBranch(fMCEvtInfoName);
  if(!fIsData && fUseGen) 
    LoadBranch(fBarrelSCName);
  if(!fIsData && fUseGen) 
    LoadBranch(fEndcapSCName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fBeamSpotName);
  LoadBranch(fPFJetName);
  LoadBranch(fTrigMaskName);
  LoadBranch(fPFMetName); 
  LoadBranch(fPhotonName);
  LoadBranch(fConversionName);
  LoadBranch(fPUEnergyDensityName);
  LoadBranch(fPFCandidateName);
  LoadBranch(fTracksName);
  if(!fIsData)
    LoadBranch(fPileupName);

  //
  // Scan generator level info for Z -> l l process 
  //
  if(fUseGen) {  
    Double_t nEcms = 0;
    UInt_t nInit = 0;
    vector<Double_t> p4tot(4,0.);
    const Double_t tolerance = 0.5;
    
    UInt_t npho = 0;
    Double_t phoPtMax = 0;
    if(fParticles) {
      const MCParticle *boson=0, *lepton1=0, *lepton2=0, *pho=0;
      for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
        const MCParticle *p = fParticles->At(i);
      
        if(nInit<3 && p->Status()==3 && p->PdgId()==2212) {
	  nInit++;
	  nEcms += p->E();
	}
	if(p->Status()==1) {
	  p4tot[0] += p->Px();
	  p4tot[1] += p->Py();
	  p4tot[2] += p->Pz();
	  p4tot[3] += sqrt((p->Px())*(p->Px()) +
	                   (p->Py())*(p->Py()) +
			   (p->Pz())*(p->Pz()) +
			   (p->Mass())*(p->Mass()));
	}
	
	if(fFSRMode==0) {
          //--------------- PYTHIA FSR mode ---------------//
          // a "branching" in the process tree is created
          // for every physical process; need to scan down
          // the lepton branches and search for photons
          //
          if( (p->PdgId() == 23) && (p->Status() == 3) ) {
            boson = p;
        
            // Loop through daughters and look for electrons      
            for(UInt_t j=0; j<boson->NDaughters(); j++) {
              const MCParticle *d = boson->Daughter(j);
              if( abs(d->PdgId())==11 || abs(d->PdgId())==13 || abs(d->PdgId())==15 ) {
            
                while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1)) {
                  const MCParticle *g = d->FindDaughter(22,kTRUE);
                  if(g && (g->Status()==1)) {
                    npho++;
                    if(g->Pt() > phoPtMax) {
                      phoPtMax = g->Pt();
                      pho = g;
                    }
                  }
                  d = d->FindDaughter(d->PdgId(),kTRUE);
                }     

                if(d->PdgId() > 0) lepton1 = d;
                if(d->PdgId() < 0) lepton2 = d;
              }  
            }                     
          }
        } else if(fFSRMode==1) {
          //--------------- HORACE FSR mode ---------------//
          // leptons and photons are daughters of the Z
          //
          if( (p->PdgId() == 23) && (p->Status() == 3) ) {
            boson = p;
          
            // Loop through daughters and look for electrons and photons
            for(UInt_t j=0; j<boson->NDaughters(); j++) {
              const MCParticle *d = boson->Daughter(j);
              if(d->PdgId()== 11 || d->PdgId()== 13 || d->PdgId()== 15) lepton1 = d;
              if(d->PdgId()==-11 || d->PdgId()==-13 || d->PdgId()==-15) lepton2 = d;
              if(d->PdgId()== 22) {
                npho++;
                if(d->Pt() > phoPtMax) {
                  phoPtMax = d->Pt();
                  pho = d;
                }
              }
            }
          }
        } else if(fFSRMode==2) {
          //--------------- PHOTOS FSR mode ---------------//
          // HACKY PHOTOS SETUP!!!!
          // Z boson has status=2
          // No daughter links (eliminated to solve SIM
          // issues), assume Z daughters are the next two
          // entries in the GEN particle list. FSR photons
          // are "found" at the end of the list of status=1
          // particles.
          //
          if( (p->PdgId() == 23) && (p->Status() == 2) ) {
            boson = p;
            const MCParticle *d = fParticles->At(i+1);
            if(d->PdgId()== 11 || d->PdgId()== 13 || d->PdgId()== 15) lepton1 = d;
            if(d->PdgId()==-11 || d->PdgId()==-13 || d->PdgId()==-15) lepton2 = d;
            d = fParticles->At(i+2);
            if(d->PdgId()== 11 || d->PdgId()== 13 || d->PdgId()== 15) lepton1 = d;
            if(d->PdgId()==-11 || d->PdgId()==-13 || d->PdgId()==-15) lepton2 = d;
        
            // a photon is an FSR photon if including it makes
            // the system closer to the Z boson mass
          
            assert(boson && lepton1 && lepton2);
            // find index of first status=-99 particle
            UInt_t ii=0;
            const MCParticle *tmpp = fParticles->At(ii);
            while(tmpp && (tmpp->Status()!=-99)) {
              ii++;
              tmpp = fParticles->At(ii);
            }
            // scan photons in "reverse" order
            const MCParticle *gamma = fParticles->At(ii-1);
            if((gamma->PdgId()==22) && (gamma->Status()==1)) {        
              FourVectorM vBoson = boson->Mom();                       // boson mass
              FourVectorM vDiEle = lepton1->Mom() + lepton2->Mom();    // dielectron mass
              Double_t dm = fabs(vDiEle.M()-vBoson.M());               // mass difference
              Double_t sm = vDiEle.M()+vBoson.M();                     // mass sum (sets order of magnitude for comparison)
              const Double_t epsilon = 0.000005;                       // defines "equals" for floating point numbers
              while(gamma && (dm > sm*epsilon)) {
                FourVectorM test = vDiEle + gamma->Mom();
                if(fabs(test.M()-vBoson.M()) < dm) {
                  vDiEle = test;
                  dm     = fabs(test.M()-vBoson.M());
                  sm     = test.M()+vBoson.M();
                  npho++;
                  if(gamma->Pt() > phoPtMax) {
                    phoPtMax = gamma->Pt();
                    pho = gamma;
                  }
              
                  gamma = fParticles->At(ii-1-npho);
                  if(gamma==0 || gamma->PdgId()!=22 || gamma->Status()!=1) 
                    gamma = 0;
              
                } else {
                  gamma = 0;
                }
              }
            }
          } 

        } else if(fFSRMode==-1) {
          //--------------- ttbar mode ---------------//
          // leptons => t quarks, boson => di-top system (A bit hacky)
          //
	  if( (p->PdgId()== 6) && (p->Status()==3) ) lepton1 = p;
	  if( (p->PdgId()==-6) && (p->Status()==3) ) lepton2 = p;
	  
	  if(!boson && lepton1 && lepton2) {
	    boson = new MCParticle(lepton1->Px()+lepton2->Px(),
	                           lepton1->Py()+lepton2->Py(),
				   lepton1->Pz()+lepton2->Pz(),
				   lepton1->E() +lepton2->E(),
				   0,0);
	  
	  }
	}
      }
      
      // Filter unphysical events
      //if(fabs(p4tot[0])>tolerance || fabs(p4tot[1]>tolerance) || fabs(p4tot[2]>tolerance) || fabs(p4tot[3] - nEcms)>tolerance)
      //  return;
      
      if(boson && lepton1 && lepton2) {         
        FillGenInfo(boson, lepton1, lepton2, pho, npho);  // fill the data structure
        fGenInfo.pxtot  = p4tot[0];
	fGenInfo.pytot  = p4tot[1];
	fGenInfo.pztot  = p4tot[2];
	fGenInfo.deltaE = p4tot[3] - nEcms;
	
	if(fFSRMode==-1) {
	  delete boson;
	  boson=0;
	}
      }
    }  
  }
  
  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  //
  ULong_t trigbits=0;
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0))
    return;
  
  IncNEventsProcessed();
  
  fMuonArr->Clear();
  fElectronArr->Clear();
  fDielectronArr->Clear();
  fPFJetArr->Clear();
  fPhotonArr->Clear();
  fPVArr->Clear();


  //
  // Get beam spot. If no beam spot information is available, default the coordinates to 99999
  //
  Double_t bsx=99999, bsy=99999, bsz=99999;
  if(fBeamSpot) {
    if(fBeamSpot->GetEntries() > 1) 
      std::cout << "********** More than 1 beam spot! **********" << std::endl;
    const BeamSpot *bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }

  //
  // Get primary vertices
  // Assumes primary vertices are ordered by sum-pT^2 (as should be in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  fVertex = 0;
  UInt_t nvertices = 0;
  Bool_t hasGoodPV = kFALSE;  
  for(UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    const Vertex *pv = fPrimVerts->At(i);
    
    // Select best PV for corrected d0; if no PV passing cuts, the first PV in the collection will be used
    //if(!pv->IsValid()) continue;
    if(pv->NTracksFit()     < fMinNTracksFit) continue;
    if(pv->Ndof()	    < fMinNdof)	      continue;
    if(fabs(pv->Z())	    > fMaxAbsZ)	      continue;
    if(pv->Position().Rho() > fMaxRho)	      continue;    
    hasGoodPV = kTRUE;
    
    nvertices++;
    FillPV(pv);
    if(!fVertex) fVertex = pv;
  }
  if(!fVertex) fVertex = fPrimVerts->At(0);
   
  //
  // Separate PF candidates into those associated
  // with PV and those associated with PU
  //
  separatePileUp(kTRUE);
   
  //
  // Loop through muons (and general tracks if desired).
  //
  vector<const Muon*> muonv;    // array of pointers to preselected muons ... 
  
  assert(fMuons);
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    
    // Use tracker tracks for kinematics when available
    const Track *muTrk=0;
    if(mu->HasTrackerTrk())     { muTrk = mu->TrackerTrk(); }
    else if(mu->HasGlobalTrk()) { muTrk = mu->GlobalTrk(); }
    else                        { muTrk = mu->StandaloneTrk(); } 
              
    if((muTrk->Eta() < fMuEtaMin) || (muTrk->Eta() > fMuEtaMax)) continue;   
    if((muTrk->Pt()  > fMuPtMax)  || (muTrk->Pt()  < fMuPtMin))  continue;
    
    FillMuon(mu);  
  }
       

  //
  // Loop through electrons.
  //
  vector<const Electron*> elev;  // array of pointers to preselected electrons ... 
  vector<double> scalev;         // array of electron energy corrections
  ElectronTools eleTools;        // helper class for electron ID decisions
  
  assert(fElectrons);                
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
   
    pair<double,double> result = fEleCorr.evaluate(ele,
                                                   fIs2011 ? fPUEnergyDensity->At(0)->RhoKt6PFJetsForIso25() : fPUEnergyDensity->At(0)->RhoKt6PFJets(),
                                                   nvertices,
                                                   GetEventHeader()->RunNum(),
                                                   false);

    double newPt = result.first*TMath::Sin(ele->Theta());
    double scale = newPt/(ele->Pt());
   
    if((ele->Pt()  < fEleEtMin)  && (newPt      < fEleEtMin))  continue;  // electron min pT cut
    if((ele->Pt()  > fEleEtMax)  && (newPt      > fEleEtMax))  continue;  // electron max pT cut
    if((ele->Eta() < fEleEtaMin) || (ele->Eta() > fEleEtaMax)) continue;  // electron eta cut
    if(!eleTools.PassSpikeRemovalFilter(ele))                  continue;  // spike cleaning
        
    elev.push_back(ele);
    scalev.push_back(scale);
    FillElectron(ele,scale);  // fill electron data object
  }
  
  // Reconstruct dielectrons
  for(UInt_t i=0; i<elev.size(); i++) {
    for(UInt_t j=i+1; j<elev.size(); j++) {
      
      // check Z mass window
      FourVectorM vDiEle = elev[i]->Mom() + elev[j]->Mom();
      if((vDiEle.M() < fMassMin) || (vDiEle.M() > fMassMax)) continue;
      
      FillDielectron(elev[i],scalev[i],elev[j],scalev[j]);  // fill dielectron data object                                    
    }
  }
  
  //
  // Loop through jets
  //
  assert(fPFJets);
  for(UInt_t i=0; i<fPFJets->GetEntries(); ++i) {
    const PFJet *jet = fPFJets->At(i);
     
    const FourVectorM rawMom = jet->RawMom();
    fJetCorrector->setJetEta(rawMom.Eta());
    fJetCorrector->setJetPt(rawMom.Pt());
    fJetCorrector->setJetPhi(rawMom.Phi());
    fJetCorrector->setJetE(rawMom.E());
    fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoKt6PFJets());
    fJetCorrector->setJetA(jet->JetArea());
    fJetCorrector->setJetEMF(-99.0);     
    
    if( ((jet->RawMom().Pt())*(fJetCorrector->getCorrection()) > fJetPtMin)
         || (jet->TrackCountingHighEffBJetTagsDisc()!=-100) ) {
       FillJet(jet);
    }
  }  
    
  //
  // Loop through photons
  //
  assert(fPhotons);
  for(UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *pho= fPhotons->At(i);
    if(pho->SCluster()->Et() > fPhotonEtMin) { FillPhoton(pho); }
  } 

/*
  //
  // Loop through superclusters
  //  
  assert(fBarrelSC);
  for(UInt_t i=0; i<fBarrelSC->GetEntries(); ++i) {
    const SuperCluster *sc = fBarrelSC->At(i);
    TVector3 scPos; scPos.SetPtEtaPhi(sc->Point().Rho(),sc->Eta(),sc->Phi());
    TVector3 vtx;   vtx.SetXYZ(fPrimVerts->At(0)->X(), fPrimVerts->At(0)->Y(), fPrimVerts->At(0)->Z());
    scPos -= vtx;
    Double_t et = sc->Energy()*(scPos.Perp())/(scPos.Mag());
    if(et > fPhotonEtMin) { FillPhoton(sc); }
  }
  assert(fEndcapSC);
  for(UInt_t i=0; i<fEndcapSC->GetEntries(); ++i) {
    const SuperCluster *sc = fEndcapSC->At(i);
    TVector3 scPos; scPos.SetPtEtaPhi(sc->Point().Rho(),sc->Eta(),sc->Phi());
    TVector3 vtx;   vtx.SetXYZ(fPrimVerts->At(0)->X(), fPrimVerts->At(0)->Y(), fPrimVerts->At(0)->Z());
    scPos -= vtx;
    Double_t et = sc->Energy()*(scPos.Perp())/(scPos.Mag());
    if(et > fPhotonEtMin) { FillPhoton(sc); }     
  } 
*/
  
  //
  // Compute MET
  //
  TLorentzVector pfmet; pfmet.SetPxPyPzE(fPFMet->At(0)->Mex(),fPFMet->At(0)->Mey(),0,0);
  Double_t trkMetx=0, trkMety=0, trkSumET=0;
  
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const Double_t trkDzCut  = 0.1;
    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*fVertex))<trkDzCut)) ||
        (pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(*fVertex))<trkDzCut)) ) {
      
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumET += pfcand->Pt();
    }
  }
  TLorentzVector trkmet;    trkmet.SetPxPyPzE(trkMetx,trkMety,0,0);

  //
  // Identify bunch crossing 0 in PU info
  //
  Int_t ibx0=-1, ibxm=-1, ibxp=-1;
  if(!fIsData) {
    assert(fPileup);    
    for(UInt_t i=0; i<fPileup->GetEntries(); ++i) {
      if(fPileup->At(i)->GetBunchCrossing() == 0) ibx0=i;
      if(fPileup->At(i)->GetBunchCrossing() ==-1) ibxm=i;
      if(fPileup->At(i)->GetBunchCrossing() == 1) ibxp=i;
    }
  }
          
  //
  // Fill event info tree
  //    
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nPU          = (ibx0>-1) ? fPileup->At(ibx0)->GetPU_NumInteractions() : 0;
  fEventInfo.nPUminus     = (ibxm>-1) ? fPileup->At(ibxm)->GetPU_NumInteractions() : 0;
  fEventInfo.nPUplus      = (ibxp>-1) ? fPileup->At(ibxp)->GetPU_NumInteractions() : 0;
  fEventInfo.nPUmean      = (ibx0>-1) ? fPileup->At(ibx0)->GetPU_NumMean() : 0;
  fEventInfo.nPUmeanminus = (ibxm>-1) ? fPileup->At(ibxm)->GetPU_NumMean() : 0;
  fEventInfo.nPUmeanplus  = (ibxp>-1) ? fPileup->At(ibxp)->GetPU_NumMean() : 0;  
  fEventInfo.triggerBits  = trigbits;
  fEventInfo.pvx          = fVertex->X();
  fEventInfo.pvy          = fVertex->Y();
  fEventInfo.pvz          = fVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.pfMET        = pfmet.Pt();
  fEventInfo.pfMETphi     = pfmet.Phi();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.trkMET       = trkmet.Pt();
  fEventInfo.trkMETphi    = trkmet.Phi();
  fEventInfo.trkSumET     = trkSumET;
  fEventInfo.rhoLowEta    = fPUEnergyDensity->At(0)->RhoKt6PFJets();
  fEventInfo.rhoHighEta   = fPUEnergyDensity->At(0)->RhoKt6PFJets();
  fEventInfo.hasGoodPV    = hasGoodPV;
  
  // Fill the tree
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillGenInfo(const MCParticle *p, const MCParticle *lepton1, const MCParticle *lepton2,
                             const MCParticle *pho, const UInt_t npho)
{
  assert(p);
  assert(lepton1);
  assert(lepton2);
  
  // Get the status=3 daughters. Not applicable for HORACE or PHOTOS modes.
  const MCParticle *vlepton1=lepton1;
  const MCParticle *vlepton2=lepton2;
  if(fFSRMode==0) {
    while(vlepton1->Status()!=3)
      vlepton1 = vlepton1->FindMother(vlepton1->PdgId(),kTRUE);
    
    while(vlepton2->Status()!=3)
      vlepton2 = vlepton2->FindMother(vlepton2->PdgId(),kTRUE);
  }
  
  fGenInfo.npho = npho;
  
  fGenInfo.id_1   = fMCEvtInfo->Id1();
  fGenInfo.id_2   = fMCEvtInfo->Id2();
  fGenInfo.x_1    = fMCEvtInfo->X1();
  fGenInfo.x_2    = fMCEvtInfo->X2();
  fGenInfo.weight = fMCEvtInfo->Weight();
  
  fGenInfo.vmass = p->Mass();
  fGenInfo.vpt   = p->Pt();
  fGenInfo.vy    = p->Rapidity();
  fGenInfo.vphi  = p->Phi();
    
  fGenInfo.lid_1  = vlepton1->PdgId();
  fGenInfo.vpt_1  = vlepton1->Pt();
  fGenInfo.veta_1 = vlepton1->Eta();
  fGenInfo.vphi_1 = vlepton1->Phi();
  
  fGenInfo.lid_2  = vlepton2->PdgId();
  fGenInfo.vpt_2  = vlepton2->Pt();  
  fGenInfo.veta_2 = vlepton2->Eta();
  fGenInfo.vphi_2 = vlepton2->Phi(); 

  FourVectorM vDiEle = lepton1->Mom() + lepton2->Mom();
  fGenInfo.mass = vDiEle.M();
  fGenInfo.pt   = vDiEle.Pt();
  fGenInfo.y    = vDiEle.Rapidity();
  fGenInfo.phi  = vDiEle.Phi();
  
  fGenInfo.pt_1  = lepton1->Pt();
  fGenInfo.eta_1 = lepton1->Eta();
  fGenInfo.phi_1 = lepton1->Phi();
  
  fGenInfo.pt_2  = lepton2->Pt();  
  fGenInfo.eta_2 = lepton2->Eta();
  fGenInfo.phi_2 = lepton2->Phi(); 
  
  fGenInfo.phopt  = (pho) ? pho->Pt()  : -999;
  fGenInfo.phoeta = (pho) ? pho->Eta() : -999;
  fGenInfo.phophi = (pho) ? pho->Phi() : -999;
  
  fGenInfo.decx = p->DecayVertex().X();
  fGenInfo.decy = p->DecayVertex().Y();
  fGenInfo.decz = p->DecayVertex().Z();
  
  fGenInfo.scEt_1  = -1;
  fGenInfo.scEta_1 = -999;
  fGenInfo.scEt_2  = -1;
  fGenInfo.scEta_2 = -999;
  fGenInfo.scMass  = -1;
  
  const Double_t matchDR = 0.2;
  
  const SuperCluster *sc1=0;
  Double_t dr1 = matchDR;  
  if(fabs(lepton1->Eta()) < 4) {    
    assert(fBarrelSC);
    for(UInt_t i=0; i<fBarrelSC->GetEntries(); ++i) {
      const SuperCluster *sc = fBarrelSC->At(i);
      Double_t dr = MathUtils::DeltaR(sc->Phi(),sc->Eta(),lepton1->Phi(),EcalEta(lepton1->Eta(),p->DecayVertex().Z(),p->DecayVertex().Rho()));
      if(dr > dr1) continue;
      sc1 = sc;
      dr1 = dr;
    }
    assert(fEndcapSC);
    for(UInt_t i=0; i<fEndcapSC->GetEntries(); ++i) {
      const SuperCluster *sc = fEndcapSC->At(i);
      Double_t dr = MathUtils::DeltaR(sc->Phi(),sc->Eta(),lepton1->Phi(),EcalEta(lepton1->Eta(),p->DecayVertex().Z(),p->DecayVertex().Rho()));
      if(dr > dr1) continue;
      sc1 = sc;
      dr1 = dr;      
    } 
  }

  const SuperCluster *sc2=0;
  Double_t dr2 = matchDR;
  if(fabs(lepton2->Eta()) < 4) {
    assert(fBarrelSC);
    for(UInt_t i=0; i<fBarrelSC->GetEntries(); ++i) {
      const SuperCluster *sc = fBarrelSC->At(i);
      Double_t dr = MathUtils::DeltaR(sc->Phi(),sc->Eta(),lepton2->Phi(),EcalEta(lepton2->Eta(),p->DecayVertex().Z(),p->DecayVertex().Rho()));
      if(dr > dr2) continue;
      sc2 = sc;
      dr2 = dr;    
    }
    assert(fEndcapSC);
    for(UInt_t i=0; i<fEndcapSC->GetEntries(); ++i) {
      const SuperCluster *sc = fEndcapSC->At(i);
      Double_t dr = MathUtils::DeltaR(sc->Phi(),sc->Eta(),lepton2->Phi(),EcalEta(lepton2->Eta(),p->DecayVertex().Z(),p->DecayVertex().Rho()));
      if(dr > dr2) continue;
      sc2 = sc;
      dr2 = dr;      
    } 
  }

  TVector3 vtx; vtx.SetXYZ(fPrimVerts->At(0)->X(), fPrimVerts->At(0)->Y(), fPrimVerts->At(0)->Z());
  TVector3 sc1Pos, sc2Pos;
  if(sc1) { 
    sc1Pos.SetPtEtaPhi(sc1->Point().Rho(),sc1->Eta(),sc1->Phi());
    sc1Pos -= vtx;
    fGenInfo.scEt_1 = sc1->Energy()*(sc1Pos.Perp())/(sc1Pos.Mag()); 
    fGenInfo.scEta_1 = sc1->Eta(); 
  }
  if(sc2) { 
    sc2Pos.SetPtEtaPhi(sc2->Point().Rho(),sc2->Eta(),sc2->Phi());
    sc2Pos -= vtx;
    fGenInfo.scEt_2 = sc2->Energy()*(sc2Pos.Perp())/(sc2Pos.Mag()); 
    fGenInfo.scEta_2 = sc2->Eta(); 
  }
  if(sc1 && sc2) {
    const Double_t m = 0.000511;
    Double_t pt1 = sqrt((sc1->Energy())*(sc1->Energy())-m*m)*(sc1Pos.Perp())/(sc1Pos.Mag());
    Double_t pt2 = sqrt((sc2->Energy())*(sc2->Energy())-m*m)*(sc2Pos.Perp())/(sc2Pos.Mag());
    FourVectorM v1(pt1,sc1Pos.Eta(),sc1Pos.Phi(),m);
    FourVectorM v2(pt2,sc2Pos.Eta(),sc2Pos.Phi(),m);
    FourVectorM diSC = v1+v2;
    fGenInfo.scMass = diSC.M();
  }  
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillMuon(const Muon *mu)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk(); }
  else if(mu->HasGlobalTrk())     { muTrk = mu->GlobalTrk(); }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); }
  assert(muTrk);                  
  
  pMuon->pt       = muTrk->Pt();
  pMuon->ptErr    = muTrk->PtErr();
  pMuon->eta      = muTrk->Eta();
  pMuon->phi      = muTrk->Phi();
  pMuon->trkIso03 = mu->IsoR03SumPt();
  pMuon->emIso03  = mu->IsoR03EmEt();
  pMuon->hadIso03 = mu->IsoR03HadEt();  
  pMuon->d0       = muTrk->D0Corrected(*fVertex);
  pMuon->dz       = muTrk->DzCorrected(*fVertex);
  pMuon->tkNchi2  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->RChi2() : 0;
  
  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  
  pMuon->trkKink    = mu->TrkKink();
  pMuon->gblKink    = mu->GlbKink();
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  //
  // NOTE:
  // It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA if
  // chamber matches for the TK and hits for the SA share DetIDs
  //   (see hypernews thread: https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html)
  //	      
  pMuon->typeBits = 0;
  if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  if(mu->IsPFMuon())         { pMuon->typeBits |= kPFMuon; }
  
  pMuon->nTkHits      = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->NHits() : 0;
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nSeg         = mu->NSegments();
  pMuon->nMatch       = mu->NMatches();
  pMuon->hltMatchBits = MatchHLT(muTrk->Eta(),muTrk->Phi());
  
  pMuon->staPt  = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Phi() : -999;
  pMuon->trkID  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->GetUniqueID() : 0;
  
  pMuon->pfPt=0;
  pMuon->pfEta=-999;
  pMuon->pfPhi=-999;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if(mu->HasTrackerTrk() && mu->TrackerTrk() == pfcand->TrackerTrk()) {
      pMuon->pfPt  = pfcand->Pt();
      pMuon->pfEta = pfcand->Eta();
      pMuon->pfPhi = pfcand->Phi();
      break;
    }
  }
  
  Double_t chIso[5], gammaIso[5], neuHadIso[5], puIso[5];
  computePFMuonIso(mu, chIso, gammaIso, neuHadIso, puIso);
  pMuon->chIso_00_01 = chIso[0];	
  pMuon->chIso_01_02 = chIso[1];
  pMuon->chIso_02_03 = chIso[2];
  pMuon->chIso_03_04 = chIso[3];
  pMuon->chIso_04_05 = chIso[4];
  pMuon->gammaIso_00_01 = gammaIso[0]; 
  pMuon->gammaIso_01_02 = gammaIso[1];
  pMuon->gammaIso_02_03 = gammaIso[2];
  pMuon->gammaIso_03_04 = gammaIso[3];
  pMuon->gammaIso_04_05 = gammaIso[4];
  pMuon->neuHadIso_00_01 = neuHadIso[0];
  pMuon->neuHadIso_01_02 = neuHadIso[1];
  pMuon->neuHadIso_02_03 = neuHadIso[2];
  pMuon->neuHadIso_03_04 = neuHadIso[3];
  pMuon->neuHadIso_04_05 = neuHadIso[4];
  pMuon->puIso_00_01 = puIso[0]; 
  pMuon->puIso_01_02 = puIso[1];
  pMuon->puIso_02_03 = puIso[2];
  pMuon->puIso_03_04 = puIso[3];
  pMuon->puIso_04_05 = puIso[4];
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillElectron(const Electron *ele, const double scale)
{
  assert(ele);
  
  TClonesArray &rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const Int_t index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron *pElectron = (TElectron*)rElectronArr[index];
                    
  pElectron->pt              = scale*(ele->Pt());
  pElectron->ptUncorr        = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->d0              = ele->BestTrk()->D0Corrected(*fVertex);
  pElectron->dz              = ele->BestTrk()->DzCorrected(*fVertex);  
  pElectron->scEt            = scale*(ele->SCluster()->Energy())*(ele->Pt())/(ele->P());
  pElectron->scEtUncorr      = (ele->SCluster()->Energy())*(ele->Pt())/(ele->P());
  pElectron->scEta           = ele->SCluster()->Eta();
  pElectron->scPhi           = ele->SCluster()->Phi();
  pElectron->r9              = ele->SCluster()->R9();
  pElectron->ecalE           = ele->EcalEnergy();
  pElectron->HoverE          = ele->HadronicOverEm();
  pElectron->EoverP          = ele->ESuperClusterOverP();
  pElectron->fBrem           = ele->FBrem();
  pElectron->deltaEtaIn      = ele->DeltaEtaSuperClusterTrackAtVtx();
  pElectron->deltaPhiIn      = ele->DeltaPhiSuperClusterTrackAtVtx();
  pElectron->sigiEtaiEta     = ele->CoviEtaiEta();
  pElectron->nExpHitsInner   = ele->BestTrk()->NExpectedHitsInner();
  pElectron->partnerDeltaCot = ele->ConvPartnerDCotTheta();
  pElectron->partnerDist     = ele->ConvPartnerDist();
  pElectron->q               = ele->Charge(); 
  
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi());
  pElectron->scID            = ele->SCluster()->GetUniqueID();
  pElectron->trkID           = (ele->HasTrackerTrk()) ? ele->TrackerTrk()->GetUniqueID() : 0;
  pElectron->isConv          = IsConversion(ele);
  
  pElectron->typeBits=0;
  if(ele->IsEcalDriven())    { pElectron->typeBits |= kEcalDriven; }
  if(ele->IsTrackerDriven()) { pElectron->typeBits |= kTrackerDriven; }
    
  pElectron->pfPt=0;
  pElectron->pfEta=-999;
  pElectron->pfPhi=-999;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFCandidates->At(i);
 
    if( (pfcand->HasTrackerTrk() && ele->TrackerTrk() == pfcand->TrackerTrk()) ||
        (pfcand->HasGsfTrk() && ele->GsfTrk() == pfcand->GsfTrk()) ) {
          pElectron->pfPt  = pfcand->Pt();
          pElectron->pfEta = pfcand->Eta();
	  pElectron->pfPhi = pfcand->Phi();
          break;
    }	 
  }
  
  pElectron->mva = fEleMVA->MVAValue(ele, fVertex);
  
  Double_t chIso[5], gammaIso[5], neuHadIso[5];
  computePFEleIso(ele, chIso, gammaIso, neuHadIso);
  pElectron->chIso_00_01 = chIso[0];	
  pElectron->chIso_01_02 = chIso[1];
  pElectron->chIso_02_03 = chIso[2];
  pElectron->chIso_03_04 = chIso[3];
  pElectron->chIso_04_05 = chIso[4];
  pElectron->gammaIso_00_01 = gammaIso[0]; 
  pElectron->gammaIso_01_02 = gammaIso[1];
  pElectron->gammaIso_02_03 = gammaIso[2];
  pElectron->gammaIso_03_04 = gammaIso[3];
  pElectron->gammaIso_04_05 = gammaIso[4];
  pElectron->neuHadIso_00_01 = neuHadIso[0];
  pElectron->neuHadIso_01_02 = neuHadIso[1];
  pElectron->neuHadIso_02_03 = neuHadIso[2];
  pElectron->neuHadIso_03_04 = neuHadIso[3];
  pElectron->neuHadIso_04_05 = neuHadIso[4];
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillDielectron(const Electron *ele1, const double scale1, const Electron *ele2, const double scale2)
{
  assert(ele1);
  assert(ele2);
  
  TClonesArray &rDielectronArr = *fDielectronArr;
  assert(rDielectronArr.GetEntries() < rDielectronArr.GetSize());
  const Int_t index = rDielectronArr.GetEntries();  
  new(rDielectronArr[index]) TDielectron();
  TDielectron *pDielectron = (TDielectron*)rDielectronArr[index];
  
  const Electron *leading = (ele1->Pt() > ele2->Pt()) ? ele1 : ele2;
  const Electron *lagging = (ele1->Pt() > ele2->Pt()) ? ele2 : ele1;
  const double scaleLeading = (ele1->Pt() > ele2->Pt()) ? scale1 : scale2;
  const double scaleLagging = (ele1->Pt() > ele2->Pt()) ? scale2 : scale1;

  FourVectorM vEle1; vEle1.SetCoordinates(scale1*(ele1->Pt()), ele1->Eta(), ele1->Phi(), ele1->Mass());
  FourVectorM vEle2; vEle2.SetCoordinates(scale2*(ele2->Pt()), ele2->Eta(), ele2->Phi(), ele2->Mass());
  FourVectorM vDiEle = vEle1+vEle2;

  pDielectron->mass = vDiEle.M(); 
  pDielectron->pt   = vDiEle.Pt(); 
  pDielectron->y    = vDiEle.Rapidity(); 
  pDielectron->phi  = vDiEle.Phi();   

  // leading electron
  pDielectron->pt_1              = scaleLeading*(leading->Pt());
  pDielectron->ptUncorr_1        = leading->Pt();
  pDielectron->eta_1             = leading->Eta();
  pDielectron->phi_1             = leading->Phi();
  pDielectron->trkIso03_1        = leading->TrackIsolationDr03();
  pDielectron->emIso03_1         = leading->EcalRecHitIsoDr03();
  pDielectron->hadIso03_1        = leading->HcalTowerSumEtDr03();
  pDielectron->d0_1              = leading->BestTrk()->D0Corrected(*fVertex);
  pDielectron->dz_1              = leading->BestTrk()->DzCorrected(*fVertex);
  pDielectron->scEt_1            = scaleLeading*(leading->SCluster()->Energy())*(leading->Pt())/(leading->P());
  pDielectron->scEtUncorr_1      = (leading->SCluster()->Energy())*(leading->Pt())/(leading->P());
  pDielectron->scEta_1           = leading->SCluster()->Eta();
  pDielectron->scPhi_1           = leading->SCluster()->Phi();
  pDielectron->r9_1              = leading->SCluster()->R9();
  pDielectron->ecalE_1           = leading->EcalEnergy();
  pDielectron->HoverE_1          = leading->HadronicOverEm();
  pDielectron->EoverP_1          = leading->ESuperClusterOverP();
  pDielectron->fBrem_1           = leading->FBrem(); 
  pDielectron->deltaEtaIn_1      = leading->DeltaEtaSuperClusterTrackAtVtx();
  pDielectron->deltaPhiIn_1      = leading->DeltaPhiSuperClusterTrackAtVtx();
  pDielectron->sigiEtaiEta_1     = leading->CoviEtaiEta();
  pDielectron->nExpHitsInner_1   = leading->BestTrk()->NExpectedHitsInner();
  pDielectron->partnerDeltaCot_1 = leading->ConvPartnerDCotTheta();
  pDielectron->partnerDist_1     = leading->ConvPartnerDist();
  pDielectron->q_1               = leading->Charge();
  pDielectron->isConv_1          = IsConversion(leading); 
  pDielectron->hltMatchBits_1    = MatchHLT(leading->SCluster()->Eta(),leading->SCluster()->Phi()); 
  pDielectron->scID_1            = leading->SCluster()->GetUniqueID();
      
  pDielectron->typeBits_1 = 0;
  if(ele1->IsEcalDriven())    { pDielectron->typeBits_1 |= kEcalDriven; }
  if(ele1->IsTrackerDriven()) { pDielectron->typeBits_1 |= kTrackerDriven; }

  // lagging electron
  pDielectron->pt_2              = scaleLagging*(lagging->Pt());
  pDielectron->ptUncorr_2        = lagging->Pt();
  pDielectron->eta_2             = lagging->Eta();
  pDielectron->phi_2             = lagging->Phi();
  pDielectron->trkIso03_2        = lagging->TrackIsolationDr03();
  pDielectron->emIso03_2         = lagging->EcalRecHitIsoDr03();
  pDielectron->hadIso03_2        = lagging->HcalTowerSumEtDr03();
  pDielectron->d0_2              = lagging->BestTrk()->D0Corrected(*fVertex);
  pDielectron->dz_2              = lagging->BestTrk()->DzCorrected(*fVertex);
  pDielectron->scEt_2            = scaleLagging*(lagging->SCluster()->Energy())*(lagging->Pt())/(lagging->P());
  pDielectron->scEtUncorr_2      = (lagging->SCluster()->Energy())*(lagging->Pt())/(lagging->P());
  pDielectron->scEta_2           = lagging->SCluster()->Eta();
  pDielectron->scPhi_2           = lagging->SCluster()->Phi();
  pDielectron->r9_2              = lagging->SCluster()->R9();
  pDielectron->ecalE_2           = lagging->EcalEnergy();
  pDielectron->HoverE_2          = lagging->HadronicOverEm();
  pDielectron->EoverP_2          = lagging->ESuperClusterOverP();
  pDielectron->fBrem_2           = lagging->FBrem(); 
  pDielectron->deltaEtaIn_2      = lagging->DeltaEtaSuperClusterTrackAtVtx();
  pDielectron->deltaPhiIn_2      = lagging->DeltaPhiSuperClusterTrackAtVtx();
  pDielectron->sigiEtaiEta_2     = lagging->CoviEtaiEta();
  pDielectron->nExpHitsInner_2   = lagging->BestTrk()->NExpectedHitsInner();
  pDielectron->partnerDeltaCot_2 = lagging->ConvPartnerDCotTheta();
  pDielectron->partnerDist_2     = lagging->ConvPartnerDist();
  pDielectron->q_2               = lagging->Charge();
  pDielectron->isConv_2          = IsConversion(lagging);
  pDielectron->hltMatchBits_2    = MatchHLT(lagging->SCluster()->Eta(),lagging->SCluster()->Phi());
  pDielectron->scID_2            = lagging->SCluster()->GetUniqueID(); 

  pDielectron->typeBits_2 = 0;
  if(ele2->IsEcalDriven())    { pDielectron->typeBits_2 |= kEcalDriven; }
  if(ele2->IsTrackerDriven()) { pDielectron->typeBits_2 |= kTrackerDriven; }

  pDielectron->pfPt_1=0;
  pDielectron->pfEta_1=-999;
  pDielectron->pfPhi_1=-999;
  pDielectron->pfPt_2=0;
  pDielectron->pfEta_2=-999;
  pDielectron->pfPhi_2=-999;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFCandidates->At(i);
 
    if( (pfcand->HasTrackerTrk() && leading->TrackerTrk() == pfcand->TrackerTrk()) ||
        (pfcand->HasGsfTrk() && leading->GsfTrk() == pfcand->GsfTrk()) ) {
          pDielectron->pfPt_1  = pfcand->Pt();
          pDielectron->pfEta_1 = pfcand->Eta();
	  pDielectron->pfPhi_1 = pfcand->Phi();
          break;
    }
    
    if( (pfcand->HasTrackerTrk() && lagging->TrackerTrk() == pfcand->TrackerTrk()) ||
        (pfcand->HasGsfTrk() && lagging->GsfTrk() == pfcand->GsfTrk()) ) {
          pDielectron->pfPt_2  = pfcand->Pt();
          pDielectron->pfEta_2 = pfcand->Eta();
	  pDielectron->pfPhi_2 = pfcand->Phi();
          break;
    }    	 
  }
  
  pDielectron->mva_1 = fEleMVA->MVAValue(leading, fVertex);
  pDielectron->mva_2 = fEleMVA->MVAValue(lagging, fVertex);
  
  Double_t chIso[5], gammaIso[5], neuHadIso[5];
  
  computePFEleIso(leading, chIso, gammaIso, neuHadIso);  
  pDielectron->chIso_00_01_1 = chIso[0];	
  pDielectron->chIso_01_02_1 = chIso[1];
  pDielectron->chIso_02_03_1 = chIso[2];
  pDielectron->chIso_03_04_1 = chIso[3];
  pDielectron->chIso_04_05_1 = chIso[4];
  pDielectron->gammaIso_00_01_1 = gammaIso[0]; 
  pDielectron->gammaIso_01_02_1 = gammaIso[1];
  pDielectron->gammaIso_02_03_1 = gammaIso[2];
  pDielectron->gammaIso_03_04_1 = gammaIso[3];
  pDielectron->gammaIso_04_05_1 = gammaIso[4];
  pDielectron->neuHadIso_00_01_1 = neuHadIso[0];
  pDielectron->neuHadIso_01_02_1 = neuHadIso[1];
  pDielectron->neuHadIso_02_03_1 = neuHadIso[2];
  pDielectron->neuHadIso_03_04_1 = neuHadIso[3];
  pDielectron->neuHadIso_04_05_1 = neuHadIso[4];
  
  computePFEleIso(lagging, chIso, gammaIso, neuHadIso);  
  pDielectron->chIso_00_01_2 = chIso[0];	
  pDielectron->chIso_01_02_2 = chIso[1];
  pDielectron->chIso_02_03_2 = chIso[2];
  pDielectron->chIso_03_04_2 = chIso[3];
  pDielectron->chIso_04_05_2 = chIso[4];
  pDielectron->gammaIso_00_01_2 = gammaIso[0]; 
  pDielectron->gammaIso_01_02_2 = gammaIso[1];
  pDielectron->gammaIso_02_03_2 = gammaIso[2];
  pDielectron->gammaIso_03_04_2 = gammaIso[3];
  pDielectron->gammaIso_04_05_2 = gammaIso[4];
  pDielectron->neuHadIso_00_01_2 = neuHadIso[0];
  pDielectron->neuHadIso_01_02_2 = neuHadIso[1];
  pDielectron->neuHadIso_02_03_2 = neuHadIso[2];
  pDielectron->neuHadIso_03_04_2 = neuHadIso[3];
  pDielectron->neuHadIso_04_05_2 = neuHadIso[4];  
}
  
//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillJet(const PFJet *jet)
{
  TClonesArray &rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const Int_t index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) TJet();
  TJet *pPFJet = (TJet*)rPFJetArr[index]; 
  
  const FourVectorM rawMom = jet->RawMom();
  fJetCorrector->setJetEta(rawMom.Eta());
  fJetCorrector->setJetPt(rawMom.Pt());
  fJetCorrector->setJetPhi(rawMom.Phi());
  fJetCorrector->setJetE(rawMom.E());
  fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoHighEta());
  fJetCorrector->setJetA(jet->JetArea());
  fJetCorrector->setJetEMF(-99.0);

  pPFJet->pt    = (rawMom.Pt())*(fJetCorrector->getCorrection());
  pPFJet->eta   = rawMom.Eta();
  pPFJet->phi   = rawMom.Phi();
  pPFJet->mass  = jet->Mass();
  pPFJet->rawPt = rawMom.Pt();
  pPFJet->area  = jet->JetArea();
  
  pPFJet->tche = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->tchp = jet->TrackCountingHighPurBJetTagsDisc();
  
  pPFJet->mcFlavor = jet->MatchedMCFlavor();
  
  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());
  
  Double_t numer=0, denom=0;
  for(UInt_t ipf=0; ipf<jet->NPFCands(); ipf++) {
    const PFCandidate *pfcand = jet->PFCand(ipf);
    if(!pfcand->HasTrk()) continue;
    numer += (pfcand->Pt())*(pfcand->Pt())*(pfcand->BestTrk()->DzCorrected(*fVertex));
    denom += (pfcand->Pt())*(pfcand->Pt());
  }
  pPFJet->dz = (denom>0) ? numer/denom : 999;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillPhoton(const Photon *pho)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pPhoton = (TPhoton*)rPhotonArr[index];
  
  pPhoton->pt		= pho->Pt(); 
  pPhoton->eta  	= pho->Eta();
  pPhoton->phi  	= pho->Phi();
  pPhoton->scEt		= pho->SCluster()->Et(); 
  pPhoton->scEta  	= pho->SCluster()->Eta();
  pPhoton->scPhi  	= pho->SCluster()->Phi();
  pPhoton->trkIso04     = pho->HollowConeTrkIsoDr04();
  pPhoton->trkIso04NoPV = IsolationTools::TrackIsolationNoPV(pho, fBeamSpot->At(0), 0.4, 0.04, 0.0, 0.015, 0.1, TrackQuality::highPurity, fTracks);
  pPhoton->emIso04      = pho->EcalRecHitIsoDr04();
  pPhoton->hadIso04	= pho->HcalTowerSumEtDr04();
  pPhoton->HoverE	= pho->HadOverEm();
  pPhoton->R9		= pho->R9();
  pPhoton->sigiEtaiEta  = pho->CoviEtaiEta();
  pPhoton->sigiPhiiPhi  = sqrt(pho->SCluster()->Seed()->CoviPhiiPhi());
  pPhoton->hltMatchBits = MatchHLT(pho->SCluster()->Eta(),pho->SCluster()->Phi());
  pPhoton->scID         = pho->SCluster()->GetUniqueID();
  pPhoton->hasPixelSeed = pho->HasPixelSeed();
  
  pPhoton->passEleVetoConvRec = PhotonTools::PassElectronVetoConvRecovery(pho,fElectrons,fConversions,fBeamSpot->At(0));
  pPhoton->passConvId         = PhotonTools::PassConversionId(pho,PhotonTools::MatchedConversion(pho,fConversions,fBeamSpot->At(0)));
  
  Bool_t passFilter=kTRUE;
  if((pho->SCluster()->Seed()->Energy() > 5.0) && 
     (pho->SCluster()->Seed()->EMax() / pho->SCluster()->Seed()->E3x3() > 0.95))
     passFilter=kFALSE;
  
  pPhoton->passSpikeFilter = passFilter;
}

void NtuplerMod::FillPhoton(const SuperCluster *sc)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pSC = (TPhoton*)rPhotonArr[index];

  TVector3 scPos; scPos.SetPtEtaPhi(sc->Point().Rho(),sc->Eta(),sc->Phi());
  TVector3 vtx;   vtx.SetXYZ(fVertex->X(), fVertex->Y(), fVertex->Z());
  scPos -= vtx;
    
  pSC->pt	    = sc->Energy()*(scPos.Perp())/(scPos.Mag()); 
  pSC->eta  	    = sc->Eta();
  pSC->phi  	    = sc->Phi();
  pSC->scEt	    = sc->Et(); 
  pSC->scEta  	    = sc->Eta();
  pSC->scPhi  	    = sc->Phi();
  pSC->trkIso04     = 0;
  pSC->trkIso04NoPV = 0;
  pSC->emIso04      = 0;
  pSC->hadIso04	    = 0;
  pSC->HoverE	    = sc->HadOverEm();
  pSC->R9	    = 0;
  pSC->sigiEtaiEta  = sqrt(sc->Seed()->CoviEtaiEta());
  pSC->sigiPhiiPhi  = sqrt(sc->Seed()->CoviPhiiPhi());
  pSC->hltMatchBits = MatchHLT(sc->Eta(),sc->Phi());
  pSC->scID         = sc->GetUniqueID();
  pSC->hasPixelSeed = kFALSE;
  
  pSC->passEleVetoConvRec = kFALSE;
  pSC->passConvId         = kFALSE;
  
  Bool_t passFilter=kTRUE;
  if((sc->Seed()->Energy() > 5.0) && 
     (sc->Seed()->EMax() / sc->Seed()->E3x3() > 0.95))
     passFilter=kFALSE;
  
  pSC->passSpikeFilter = passFilter;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::FillPV(const Vertex *pv) 
{
  TClonesArray &rPVArr = *fPVArr;
  assert(rPVArr.GetEntries() < rPVArr.GetSize());
  const Int_t index = rPVArr.GetEntries();  
  new(rPVArr[index]) TVertex();
  TVertex *pVertex = (TVertex*)rPVArr[index];
  
  pVertex->nTracksFit = pv->NTracksFit();
  pVertex->ndof       = pv->Ndof();      
  pVertex->chi2       = pv->Chi2();  
  pVertex->x          = pv->X();
  pVertex->y          = pv->Y();
  pVertex->z          = pv->Z();
  
  pVertex->sumPt=0;
  for(UInt_t itrk=0; itrk<pv->NTracks(); itrk++)
    pVertex->sumPt += pv->Trk(itrk)->Pt();				   
}
      
//--------------------------------------------------------------------------------------------------
ULong_t NtuplerMod::MatchHLT(const Double_t eta, const Double_t phi)
{
  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
   
      while(to) {             
        if(to->IsHLT()) {          
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        } 
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

ULong_t NtuplerMod::MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi)
{
  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  const Double_t hltMatchPtFrac = 1;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        if(to->IsHLT()) {
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

//--------------------------------------------------------------------------------------------------
Bool_t NtuplerMod::IsConversion(const Electron *ele) 
{
  Bool_t isGoodConversion = kFALSE;
  
  const UInt_t   nWrongHitsMax = 0;
  const Double_t probMin       = 1e-6;
  const Double_t lxyMin        = 2.0;
  const Bool_t   matchCkf      = kTRUE;
  const Bool_t   requireArbitratedMerged = kFALSE;
  
  for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (fConversions->At(ifc)->Daughter(d))->Trk();
      if(ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk)) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (fConversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || fConversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (fConversions->At(ifc)->LxyCorrected(fVertex) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (fConversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            const StableData *sd = dynamic_cast<const StableData*>
              (fConversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }

    if(isGoodConversion == kTRUE) break;
    
  } // loop over all fConversions 

  return isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::separatePileUp(const Bool_t checkClosestZVertex)
{
  assert(fPFPileUp);
  assert(fPFNoPileUp);

  fPFPileUp->Reset();
  fPFNoPileUp->Reset();

  for(UInt_t i=0; i<fPFCandidates->GetEntries(); i++) {
    const PFCandidate *pf = fPFCandidates->At(i);
    assert(pf);

    if(pf->PFType() == PFCandidate::eHadron) {
      if(pf->HasTrackerTrk() &&
         fVertex->HasTrack(pf->TrackerTrk()) &&
         fVertex->TrackWeight(pf->TrackerTrk()) > 0)
      {
        fPFNoPileUp->Add(pf);
      
      } else {
        Bool_t vertexFound = kFALSE;
        const Vertex *closestVtx = 0;
        Double_t dzmin = 10000;

        for(UInt_t j=0; j<fPrimVerts->GetEntries(); j++) {
          const Vertex *vtx = fPrimVerts->At(j);
          assert(vtx);

          if(pf->HasTrackerTrk() &&
             vtx->HasTrack(pf->TrackerTrk()) &&
             vtx->TrackWeight(pf->TrackerTrk()) > 0)
          {
            vertexFound = kTRUE;
            closestVtx = vtx;
            break;
          }

          Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
          if(dz < dzmin) {
            closestVtx = vtx;
            dzmin = dz;
          }
        }

        if(checkClosestZVertex) {
          // Fallback: if track is not associated with any vertex,
          // associate it with the vertex closest in z
          if(vertexFound || closestVtx != fVertex)
            fPFPileUp->Add(pf);
          else
            fPFNoPileUp->Add(pf);
        
	} else {
          if(vertexFound && closestVtx != fVertex)
            fPFPileUp->Add(pf);
          else
            fPFNoPileUp->Add(pf); // Ridiculous but that's how it is
        }
      }
    } else {
      fPFNoPileUp->Add(pf);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void NtuplerMod::computePFMuonIso(const Muon *muon, Double_t *chIso, Double_t *gammaIso, Double_t *neuHadIso, Double_t *puIso)
{
  //**********************************************************
  //Isolation variables
  //**********************************************************
  chIso[0] = 0;
  chIso[1] = 0;
  chIso[2] = 0;
  chIso[3] = 0;
  chIso[4] = 0;
  gammaIso[0] = 0;
  gammaIso[1] = 0;
  gammaIso[2] = 0;
  gammaIso[3] = 0;
  gammaIso[4] = 0;
  neuHadIso[0] = 0;
  neuHadIso[1] = 0;
  neuHadIso[2] = 0;
  neuHadIso[3] = 0;
  neuHadIso[4] = 0;
  puIso[0] = 0;
  puIso[1] = 0;
  puIso[2] = 0;
  puIso[3] = 0;
  puIso[4] = 0;
  
  // Use tracker track when available
  const Track *muTrk=0;
  if     (muon->HasTrackerTrk())    { muTrk = muon->TrackerTrk(); }
  else if(muon->HasGlobalTrk())     { muTrk = muon->GlobalTrk(); }
  else if(muon->HasStandaloneTrk()) { muTrk = muon->StandaloneTrk(); }  
  
  for(UInt_t i=0; i<fPFNoPileUp->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFNoPileUp->At(i);
    
    Double_t dr = MathUtils::DeltaR(pfcand->Phi(),pfcand->Eta(),muTrk->Phi(),muTrk->Eta());
    
    if(dr<0.5) {      
      
      if(pfcand->PFType()==PFCandidate::EPFType::eHadron) { // Charged
	
	if     (dr < 0.1) chIso[0] += pfcand->Pt();
	else if(dr < 0.2) chIso[1] += pfcand->Pt();
	else if(dr < 0.3) chIso[2] += pfcand->Pt();
	else if(dr < 0.4) chIso[3] += pfcand->Pt();
	else if(dr < 0.5) chIso[4] += pfcand->Pt();
	
      } else if(pfcand->PFType()==PFCandidate::EPFType::eGamma) { // Gamma
	if(pfcand->Pt()<0.5) continue;
	
	if     (dr < 0.1) gammaIso[0] += pfcand->Pt();
	else if(dr < 0.2) gammaIso[1] += pfcand->Pt();
	else if(dr < 0.3) gammaIso[2] += pfcand->Pt();
	else if(dr < 0.4) gammaIso[3] += pfcand->Pt();
	else if(dr < 0.5) gammaIso[4] += pfcand->Pt();
      
      } else if(pfcand->PFType()==PFCandidate::EPFType::eNeutralHadron) { // Neutral Hadron
	if(pfcand->Pt()<0.5) continue;
	
	if     (dr < 0.1) neuHadIso[0] += pfcand->Pt();
	else if(dr < 0.2) neuHadIso[1] += pfcand->Pt();
	else if(dr < 0.3) neuHadIso[2] += pfcand->Pt();
	else if(dr < 0.4) neuHadIso[3] += pfcand->Pt();
	else if(dr < 0.5) neuHadIso[4] += pfcand->Pt();
      }
    }
  }
  
  for(UInt_t i=0; i<fPFPileUp->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFPileUp->At(i);
    if(!pfcand->TrackerTrk()) continue;
    if(pfcand->Pt()<0.5)      continue;
    
    Double_t dr = MathUtils::DeltaR(pfcand->Phi(),pfcand->Eta(),muTrk->Phi(),muTrk->Eta());    
    
    if(dr < 0.5) {
      if     (dr < 0.1) puIso[0] += pfcand->Pt();
      else if(dr < 0.2) puIso[1] += pfcand->Pt();
      else if(dr < 0.3) puIso[2] += pfcand->Pt();
      else if(dr < 0.4) puIso[3] += pfcand->Pt();
      else if(dr < 0.5) puIso[4] += pfcand->Pt();
    }
  }
}

void NtuplerMod::computePFEleIso(const Electron *electron, Double_t *chIso, Double_t *gammaIso, Double_t *neuHadIso)
{
  //**********************************************************
  //Isolation variables
  //**********************************************************
  chIso[0] = 0;
  chIso[1] = 0;
  chIso[2] = 0;
  chIso[3] = 0;
  chIso[4] = 0;
  gammaIso[0] = 0;
  gammaIso[1] = 0;
  gammaIso[2] = 0;
  gammaIso[3] = 0;
  gammaIso[4] = 0;
  neuHadIso[0] = 0;
  neuHadIso[1] = 0;
  neuHadIso[2] = 0;
  neuHadIso[3] = 0;
  neuHadIso[4] = 0;
  
  for(UInt_t i=0; i<fPFNoPileUp->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFNoPileUp->At(i);
    
    if(pfcand->PFType() == PFCandidate::eGamma && 
       pfcand->MvaGamma()>0.99 && electron->CorrectedNExpectedHitsInner()>0 &&
       pfcand->SCluster() && electron->SCluster() && pfcand->SCluster() == electron->SCluster()) 
      continue;
    
    Double_t dr = MathUtils::DeltaR(pfcand->Phi(),pfcand->Eta(),electron->Phi(),electron->Eta());
    if(dr<0.5) {
      
      if(pfcand->PFType() == PFCandidate::eHadron) { // Charged hadron        
	if(electron->IsEE() && dr<0.015) continue;  // footprint removal
	
	if     (dr < 0.1) chIso[0] += pfcand->Pt();
	else if(dr < 0.2) chIso[1] += pfcand->Pt();
	else if(dr < 0.3) chIso[2] += pfcand->Pt();
	else if(dr < 0.4) chIso[3] += pfcand->Pt();
	else if(dr < 0.5) chIso[4] += pfcand->Pt();
	
      } else if(pfcand->PFType()==PFCandidate::EPFType::eGamma) { // Gamma
	if(electron->IsEE() && dr<0.08) continue;  // footprint removal
	
	if     (dr < 0.1) gammaIso[0] += pfcand->Pt();
	else if(dr < 0.2) gammaIso[1] += pfcand->Pt();
	else if(dr < 0.3) gammaIso[2] += pfcand->Pt();
	else if(dr < 0.4) gammaIso[3] += pfcand->Pt();
	else if(dr < 0.5) gammaIso[4] += pfcand->Pt();
      
      } else if(pfcand->PFType() == PFCandidate::eNeutralHadron) { // Neutral Hadron
	
	if     (dr < 0.1) neuHadIso[0] += pfcand->Pt();
	else if(dr < 0.2) neuHadIso[1] += pfcand->Pt();
	else if(dr < 0.3) neuHadIso[2] += pfcand->Pt();
	else if(dr < 0.4) neuHadIso[3] += pfcand->Pt();
	else if(dr < 0.5) neuHadIso[4] += pfcand->Pt();
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
Double_t NtuplerMod::EcalEta(const Double_t eta, const Double_t z, const Double_t rho) {
  // From https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideEcalRecoClustering
  
  const Double_t R_ECAL            = 136.5;
  const Double_t Z_ENDCAP          = 328.0;
  const Double_t ETA_BARREL_ENDCAP = 1.479;
  const Double_t pi                = 3.14159265358979;

  if(eta!= 0.) {
    Double_t theta = 0.0;
    Double_t zecal = (R_ECAL - rho)*sinh(eta) + z;
      
    if(zecal != 0.0) theta = atan(R_ECAL/zecal);
    if(theta<0.0) theta = theta + pi;

    Double_t detEta = -log(tan(0.5*theta));
      
    if(fabs(detEta) > ETA_BARREL_ENDCAP) {
      Double_t zend = Z_ENDCAP;
      if(eta<0.0) zend = -zend;
      Double_t zlen = zend - z;
      Double_t rr = zlen/sinh(eta);
      theta = atan((rr+rho)/zend);
      if(theta<0.0) theta = theta + pi;
      detEta = - log(tan(0.5*theta));
    }
    return detEta;
      
  } else {
    return eta;
  }
}
