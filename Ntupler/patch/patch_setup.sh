if test -z $CMSSW_VERSION; then
    echo "Need CMSSW project area setup!";
    exit 1;
fi

echo "Patching MitPhysics..."
cp MitPhysics/Utils/interface/ElectronEnergyRegression.h      $CMSSW_BASE/src/MitPhysics/Utils/interface/
cp MitPhysics/Utils/interface/ElectronEnergySmearingScaling.h $CMSSW_BASE/src/MitPhysics/Utils/interface/
cp MitPhysics/Utils/interface/ElectronEpCombination.h         $CMSSW_BASE/src/MitPhysics/Utils/interface/
cp MitPhysics/Utils/interface/ElectronLinearityCorrection.h   $CMSSW_BASE/src/MitPhysics/Utils/interface/
cp MitPhysics/Utils/interface/ElectronMomentumCorrection.h    $CMSSW_BASE/src/MitPhysics/Utils/interface/

cp MitPhysics/Utils/src/ElectronEnergyRegression.cc      $CMSSW_BASE/src/MitPhysics/Utils/src/
cp MitPhysics/Utils/src/ElectronEnergySmearingScaling.cc $CMSSW_BASE/src/MitPhysics/Utils/src/
cp MitPhysics/Utils/src/ElectronEpCombination.cc         $CMSSW_BASE/src/MitPhysics/Utils/src/
cp MitPhysics/Utils/src/ElectronLinearityCorrection.cc   $CMSSW_BASE/src/MitPhysics/Utils/src/
cp MitPhysics/Utils/src/ElectronMomentumCorrection.cc    $CMSSW_BASE/src/MitPhysics/Utils/src/

cp MitPhysics/data/eleEnergyRegWeights_WithSubClusters_VApr15.root $CMSSW_BASE/src/MitPhysics/data/
cp MitPhysics/data/scalesCorr.csv                                  $CMSSW_BASE/src/MitPhysics/data/
cp MitPhysics/data/smearsCorrType1.csv                             $CMSSW_BASE/src/MitPhysics/data/
cp MitPhysics/data/smearsCorrType2.csv                             $CMSSW_BASE/src/MitPhysics/data/
cp MitPhysics/data/smearsCorrType3.csv                             $CMSSW_BASE/src/MitPhysics/data/
cp MitPhysics/data/linearityNewReg-May2013.csv                     $CMSSW_BASE/src/MitPhysics/data/

cp MitPhysics/Utils/BuildFile.xml                $CMSSW_BASE/src/MitPhysics/Utils/BuildFile.xml
cp MitPhysics/Utils/src/MitPhysicsUtilsLinkDef.h $CMSSW_BASE/src/MitPhysics/Utils/src/MitPhysicsUtilsLinkDef.h
