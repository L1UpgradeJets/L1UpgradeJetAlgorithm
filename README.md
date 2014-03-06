L1UpgradeJetAlgorithm
=====================

This is the developer release for the L1 upgrade jet algorithm.

Recommended release: <b>CMSSW_6_2_2</b>

Installation instructions
-------------------------

To checkout the latest <i>development</i> release:

<code>
  cmsrel CMSSW_6_2_2<br>
  cd CMSSW_6_2_2/src<br>
  cmsenv<br>
  git clone -b developer-release-v01 git://github.com/MarkBaber/L1UpgradeJetAlgorithm.git .<br>
  scram b -j8<br>
</code>


Repositry contents
------------------

The repositry contains the following directories:

    <b>AnalyseUpgradeJets/</b>        - Contains jet analyser code
    <b>SLHCUpgradeSimulations/</b>    - Contains jet emulator code
    <b>SimDataFormats/</b>            - Contains jet emulator dataformats

The core of the emulator code is run from the config files in <i>SLHCUpgradeSimulations/L1CaloTrigger/python/</i>:
        <i>SLHCCaloTrigger_cfi.py</i>
        <i>SLHCCaloTrigger_cfi.py</i>
    
Which in turn run modules in:
        <i>SLHCUpgradeSimulations/L1CaloTrigger/plugins</i>



Running and analysing the L1 upgrade jet algorithm
--------------------------------------------------

To run the jet emulator:
    
    - Run the config file <i>ProduceJetCollections.py</i> in the <i>AnalyseUpgradeJets</i> directory:
    <code>    cmsRun ProduceJetCollections.py </code>
    - This will by default produce a ROOT file <i>JetCollections.root</i> which contains a large amount of information from: caloTowers and trigger towers to offline and online jet collections and energy sums.




