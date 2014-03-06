

#include "AnalyseUpgradeJets/AnalyzeJets/interface/histBooker.h"


void 
histBooker::book1D( TString histName, TFileDirectory f, TString histTitle, histParam xParam )
{

  // ****************************************
  // *              Validation              *
  // ****************************************

  // Check that the histogram has not already been booked in the given sub directory
  TString subDir = f.fullPath();
  std::list< TString >::iterator histIt = std::find( hist1DList[subDir].begin(), hist1DList[subDir].end(), histName );

  if ( histIt != hist1DList[subDir].end() ){ // Histogram of same name already exists
    edm::LogWarning("HistogramConflict") << "1D Histogram with name: '" << histName << "' in directory: " << subDir << " booked multiple times\n";
  }
  else{
    hist1DList[subDir].push_back( histName );
  }
  
  // Book the histogram
  (*histogram1D)[ histName ] = f.make<TH1D>( histName, histTitle, xParam.bins, xParam.low, xParam.high ); 

}

//Book a PT, Eta and Phi distribution with default parameters
void
histBooker::book1DPtEtaPhi( TString histName, TFileDirectory f, TString histTitle )
{

  book1D( histName + "_PT",  f, histTitle + " Jet p_{T};p_{T} (GeV);Entries", pOffPT);
  book1D( histName + "_Eta", f, histTitle + " Jet #eta;#eta;Entries", pOffEta);
  book1D( histName + "_Phi", f, histTitle + " Jet #phi;#phi;Entries", pOffPhi);

 
} 



void 
histBooker::bookTProfile( TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam )
{

  // ****************************************
  // *              Validation              *
  // ****************************************

  // Check that the histogram has not already been booked in the given sub directory
  TString subDir = f.fullPath();
  std::list< TString >::iterator histIt = std::find( hist1DList[subDir].begin(), hist1DList[subDir].end(), histName );

  if ( histIt != hist1DList[subDir].end() ){ // Histogram of same name already exists
    edm::LogWarning("HistogramConflict") << "1D Histogram with name: '" << histName << "' in directory: " << subDir << " booked multiple times\n";
  }
  else{
    hist1DList[subDir].push_back( histName );
  }
  
  // Book the histogram - option "s" sets error of points as RMS
  (*histogram1D)[ histName ] = f.make<TProfile>( histName, histTitle, xParam.bins, xParam.low, xParam.high, yParam.low, yParam.high, "s" ); 

}


// ********************************************************************************
// *                              Book 2D histograms                              *
// ******************************************************************************** 
void 
histBooker::book2D( TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam )
{

  // ****************************************
  // *              Validation              *
  // ****************************************

  // Check that the histogram has not already been booked in the given sub directory
  TString subDir = f.fullPath();
  std::list< TString >::iterator histIt = std::find( hist2DList[subDir].begin(), hist2DList[subDir].end(), histName );

  if ( histIt != hist2DList[subDir].end() ){ // Histogram of same name already exists
    edm::LogWarning("HistogramConflict") << "2D Histogram with name: '" << histName << "' in directory: " << subDir << " booked multiple times\n";
  }
  else{
    hist2DList[subDir].push_back( histName );
  }
  
  // Book the histogram
  (*histogram2D)[ histName ] =  f.make<TH2D>( histName, histTitle, xParam.bins, xParam.low, xParam.high, yParam.bins, yParam.low, yParam.high ); 

}



void 
histBooker::book2DTProf( TString histName, TFileDirectory f, TString histTitle, histParam xParam, histParam yParam )
{

  book2D( histName, f, histTitle, xParam, yParam );
  bookTProfile( histName + "_prof", f, histTitle, xParam, yParam );

}



