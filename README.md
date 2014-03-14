L1UpgradeJetAlgorithm
=====================

Repo for the development/calibration of L1 Upgrade Jets

Running instructions:

  -Checkout the package from git: https://github.com/L1UpgradeJets/L1UpgradeJetAlgorithm/tree/pu140

  -git clone -b pu140 git@github.com:L1UpgradeJets/L1UpgradeJetAlgorithm.git .

  -Package is AnalyseUpgradeJets

  -Run ProduceJets.py over suitable sample to produce ntuples containing jets and variables

  -Run on these with AnalyseJets_GEN.py to produce histograms for calibration

  -Edit filename in macro getCalibrationJETMET.cpp to take AnalyseJets_GEN.py output

  -Load the macro in ROOT and run getCalibration()

  -Load the macro getReCalibJETMET.C

  -Run calibrateFile(“input”,minPt,maxPt) with the output root file of getCalibration() as the input

  -Save the LUT output to the terminal

  -Put the LUT into AnalyseJets_GEN.py and run again for testing
