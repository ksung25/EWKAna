// $Id: MitPhysicsUtilsLinkDef.h,v 1.28 2012/12/07 14:12:47 mhchan Exp $

#ifndef MITPHYSICS_UTILS_LINKDEF_H
#define MITPHYSICS_UTILS_LINKDEF_H

#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MatchingTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/RecoilTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetLeptonTools.h"
#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/VertexMVA.h"
#include "MitPhysics/Utils/interface/PUReweighting.h"
#include "MitPhysics/Utils/interface/PUReweightingMulti.h"
#include "MitPhysics/Utils/interface/PhotonFix.h"
#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitPhysics/Utils/interface/TauIsoMVA.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitPhysics/Utils/interface/PFMetCorrectionTools.h"
#include "MitPhysics/Utils/interface/MVAVBF.h"
#include "MitPhysics/Utils/interface/ElectronEnergyRegression.h"
#include "MitPhysics/Utils/interface/ElectronEnergySmearingScaling.h"
#include "MitPhysics/Utils/interface/ElectronEpCombination.h"
#include "MitPhysics/Utils/interface/ElectronLinearityCorrection.h"
#include "MitPhysics/Utils/interface/ElectronMomentumCorrection.h"
#include "MitPhysics/Utils/interface/AntiElectronIDMVA3.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;


#pragma link C++ class mithep::DiTauSystem;
#pragma link C++ class mithep::GeneratorTools;
#pragma link C++ class mithep::IsolationTools;
#pragma link C++ class mithep::MatchingTools;
#pragma link C++ class mithep::MuonTools;
#pragma link C++ enum mithep::MuonTools::EMuonEffectiveAreaType;
#pragma link C++ enum mithep::MuonTools::EMuonEffectiveAreaTarget; 
#pragma link C++ class mithep::MetTools;
#pragma link C++ class mithep::ElectronTools;
#pragma link C++ class mithep::PhotonTools;
#pragma link C++ class mithep::JetTools;
#pragma link C++ class mithep::MetLeptonTools;
#pragma link C++ class mithep::MetTools;
#pragma link C++ class mithep::RecoilTools;
#pragma link C++ enum mithep::ElectronTools::EElIdType;
#pragma link C++ enum mithep::ElectronTools::EElIsoType;
#pragma link C++ enum mithep::ElectronTools::EElectronEffectiveAreaType;
#pragma link C++ enum mithep::ElectronTools::EElectronEffectiveAreaTarget;
#pragma link C++ class mithep::VertexTools;
#pragma link C++ class mithep::VertexMVA;
#pragma link C++ class mithep::PUReweighting;
#pragma link C++ class mithep::PUReweightingMulti;
#pragma link C++ class PhotonFix;
#pragma link C++ class mithep::EGEnergyCorrector; 
#pragma link C++ class mithep::ElectronIDMVA; 
#pragma link C++ class mithep::MVATools; 
#pragma link C++ class mithep::MuonIDMVA; 
#pragma link C++ class mithep::JetIDMVA; 
#pragma link C++ class mithep::TauIsoMVA; 
#pragma link C++ class mithep::MVAMet; 
#pragma link C++ class mithep::PFMetCorrectionTools;
#pragma link C++ class mithep::MVAVBF;
#pragma link C++ class mithep::ElectronEnergyRegression;
#pragma link C++ class mithep::ElectronEnergySmearingScaling;
#pragma link C++ class mithep::ElectronEpCombination;
#pragma link C++ class mithep::ElectronLinearityCorrection;
#pragma link C++ class mithep::ElectronMomentumCorrection;
#pragma link C++ class mithep::AntiElectronIDMVA3;
#endif
