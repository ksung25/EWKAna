#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TDielectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TVertex.hh"
#endif

// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimNtuples(const TString input = "skim.input") 
{
  gBenchmark->Start("SkimNtuples");
  
  TString outfilename;          // output of skimming 
  vector<TString> infilenames;  // list input ntuple files to be skimmed
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input.Data()); 
  assert(ifs.is_open());
  string line;
  getline(ifs,line); 
  outfilename = line;
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();

  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TDielectron::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();
  mithep::TVertex::Class()->IgnoreTObjectStreamer();
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr   = new TClonesArray("mithep::TElectron");
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  TClonesArray *muonArr       = new TClonesArray("mithep::TMuon");
  TClonesArray *pfJetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr     = new TClonesArray("mithep::TPhoton");
  TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");
  
  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  
  TFile* outfile = new TFile(outfilename, "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",       &info);
  outEventTree->Branch("Electron",   &electronArr);
  outEventTree->Branch("Dielectron", &dielectronArr);
  outEventTree->Branch("Muon",       &muonArr);
  outEventTree->Branch("PFJet",      &pfJetArr);
  outEventTree->Branch("Photon",     &photonArr);
  outEventTree->Branch("PV",         &pvArr);

  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Skimming " << infilenames[ifile] << "..." << endl;
    TFile *infile = new TFile(infilenames[ifile]);
    assert(infile);
    
    TTree *eventTree = (TTree*)infile->Get("Events");
    assert(eventTree);
    
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron",   &electronArr);   TBranch *electronBr   = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Dielectron", &dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("Muon",       &muonArr);       TBranch *muonBr       = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("PFJet",      &pfJetArr);      TBranch *pfJetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("Photon",     &photonArr);     TBranch *photonBr     = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PV",         &pvArr);         TBranch *pvBr         = eventTree->GetBranch("PV");
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
      infoBr->GetEntry(ientry);

      electronArr->Clear();   electronBr->GetEntry(ientry);
      dielectronArr->Clear(); dielectronBr->GetEntry(ientry);
      muonArr->Clear();       muonBr->GetEntry(ientry);
      pfJetArr->Clear();      pfJetBr->GetEntry(ientry);
      photonArr->Clear();     photonBr->GetEntry(ientry);
      pvArr->Clear();         pvBr->GetEntry(ientry);
      
      nInputEvts++;
            
      Bool_t keep = kFALSE;

      if(dielectronArr->GetEntriesFast() > 0)
        keep = kTRUE;
  
      if(keep) {
        outEventTree->Fill();
        nPassEvts++;
      }
    }
  }
  
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete electronArr;
  delete dielectronArr;
  delete muonArr;
  delete pfJetArr;
  delete photonArr;
  delete pvArr;
    
  std::cout << outfilename << " created!" << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;
  
  gBenchmark->Show("SkimNtuples");
}  
