#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include "TGraphErrors.h"
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>
#include <stdlib.h>

#include <map>
#include <set>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


// #include <string>
// #include <sys/types.h>
// #include <sys/stat.h>
/*

TODO: MAKE THIS A LIST OF VECTORS TO RUN OVER
WILL ENABLE PU-PU COMPARISONS!!! AND REDUCE WAITING TIMES BETWEEN PROCESSING PLOTS

*/


// ****************************************************************************************************
// <JetMET stuff>

/// default fit with gaussian in niter iteration of mean 
void fit_gaussian(TH1D*& hrsp,
                  const double nsigma,
                  const double fitMin,
                  const int niter);

void adjust_fitrange(TH1* h,double& min,double& max);

template <class T> 
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&)) 
{ 
  std::istringstream iss(s); 
  return !(iss >> f >> t).fail();
}



// </JetMET stuff>
// ****************************************************************************************************


//

//#define TPROFILE
//#define ANALYSE_TH1
#define ANALYSE_TH2

// Switch on whether to fit the L1 pt vs inverse response curve
//#define FITTING




double CrystalBall(double* x, double* par){
  //http://en.wikipedia.org/wiki/Crystal_Ball_function
  double xcur = x[0];
  double alpha = par[0];
  double n = par[1];
  double mu = par[2];
  double sigma = par[3];
  double N = par[4];
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
  double A; double B;
  if (alpha < 0){
    A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2);
    B = n/(-1*alpha) + alpha;}
  else {
    A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2);
    B = n/alpha - alpha;}
  double f;
  if ((xcur-mu)/sigma > (-1)*alpha)
    f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/
		    (2*sigma*sigma));
  else
    f = N*A*pow((B- (xcur-mu)/sigma),(-1*n));
  delete exp;
  return f;
}




// Structure for containing calibration data from fits
struct calibData{
  
  calibData(int _iEta, std::vector<double> _p):iEta(_iEta), p(_p){}

  int iEta;                // iEta bin
  std::vector<double> p;   // polynomial value, index = order
  
};


// Print debugging messages                                                                                                                                    
//#define VERBOSE 

// Print histograms
#define PRINT_HIST
 

#ifdef VERBOSE
#    define PRINT(outputStr) std::cout << "\n************************************************************\n" << "Making: " << (outputStr) << "\n" << "************************************************************\n\n";
#    define SUBPRINT(outputStr) std::cout << "\t" << outputStr << "\n";
#    define QUIT(outputStr) std::cout << "\n\n\n" << outputStr << "\n\n\n"; exit(0);
#else
#    define PRINT(outputStr)
#    define SUBPRINT(outputStr)
#    define QUIT(outputStr)
#endif


using namespace std;


TString calibLUTFile = "etaPtCalibrationLUT.txt";
TString calibFitFile = "etaPtCalibrationFit.txt";

// Minimum and maximum jet pT (GeV) fit range
double xMin   = 25;
double xMax   = 300;
// The OffPT_vs_L1PT fitting function
TF1 *fitOffPT_vs_L1PT = new TF1("fitOffPT_vs_L1PT", "pol2", xMin, xMax);
// Extract the number of fit parameters


// pT step in scanning fit to jet pT
double ptStep = 5;



//======================
// Output options
//=======================


// Mk2
//##############################

//TString filename = "/home/hep/mb1512/JetAnalyser/CMSSW_6_0_1_PostLS1v2_patch4/src/SingleMu_05Oct_Mk2_SMALL.root";
//TString filename = "/vols/cms04/mb1512/Batch/2013-10-24_SingleMu_05Oct/SingleMu_05Oct.root";
//TString filename = "/vols/cms04/mb1512/Batch/2013-10-28_SingleMu_18Oct_9x9/SingleMu_18Oct_9x9.root";
// TString filename = "/vols/cms04/mb1512/Batch/2013-10-28_SingleMu_18Oct_11x11/SingleMu_18Oct_11x11.root";
// TString filename = "/vols/cms04/mb1512/Batch/2013-10-28_SingleMu_18Oct_13x13/SingleMu_18Oct_13x13.root";
// TString filename = "/vols/cms04/mb1512/Batch/2013-10-28_SingleMu_18Oct_15x15/SingleMu_18Oct_15x15.root";

// NEW FILES
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-04_SingleMu_18Oct_11x11/SingleMu_18Oct_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-08_SingleMu_18Oct_11x11_rev_4/SingleMu_18Oct_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-14_SingleMu_18Oct_11x11/SingleMu_18Oct_11x11.root";

// dR = 0.3, lead jet matching
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-19_SingleMu_18Oct_11x11_rev_1/SingleMu_18Oct_11x11.root";
// dR = 0.3, 2 jet matching
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-19_SingleMu_18Oct_11x11_rev_2/SingleMu_18Oct_11x11.root";
// dR = 0.3, all jet matching
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-19_SingleMu_18Oct_11x11_rev_6/SingleMu_18Oct_11x11.root";

// dR = 0.3, NVTX binned
//TString filename = "/vols/cms04/mb1512/Batch/2013-11-20_SingleMu_18Oct_11x11_rev_1/SingleMu_18Oct_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2014-01-24_SingleMu_12Dec_11x11/SingleMu_12Dec_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2014-01-25_SingleMu_12Dec_11x11_rev_1/SingleMu_12Dec_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-03_SingleMu_12Dec_11x11_rev_1/SingleMu_12Dec_11x11.root";
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-06_SingleMu_12Dec_11x11_rev_4/SingleMu_12Dec_11x11.root";

// Full-TT granularity
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-12_SingleMu_12Dec_11x11/SingleMu_12Dec_11x11.root";

// 5 GeV GCT seed region
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-13_SingleMu_12Dec_11x11_rev_1/SingleMu_12Dec_11x11.root";
// pT -> pT/2 fix
//
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-17_SingleMu_12Dec_11x11/SingleMu_12Dec_11x11.root";
// Uncalib GCT jets
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-23_SingleMu_20Feb_11x11_rev_3/SingleMu_20Feb_11x11.root";
// Threshold fixed jets
//TString filename = "/vols/cms04/mb1512/Batch/2014-02-27_SingleMu_26Feb_11x11/SingleMu_26Feb_11x11.root";

//TString filename = "/vols/cms04/mb1512/Batch/2014-03-03_SingleMu_26Feb_11x11/SingleMu_26Feb_11x11.root";


// Phase 2
TString filename = "/vols/cms02/ace09/ntuples/trigger/jets/upgd13d/OutputJetsQcdHighPtEta2p5.root";


// Directory inside ROOT file
//   TString ROOTdir  = "analyzer/";
TString ROOTdir  = "";
// Directory to store plots
//  TString plotDirectory = "plots/Mk1Release/Presentation/8x8_PreFix/";
TString plotDirectory = "/vols/cms02/ace09/trigger/jets/upgrade_140pu/CMSSW_6_2_2/src/AnalyseUpgradeJets/CalibrateQcd/final";



//========================
// Function prototypes
//========================

std::vector <TString> splitString(TString inputStr, TString splitChar);
Bool_t objectExists(TFile* file,  TString directory, TString objectName);


// Returns a vector [ntuplePath, objectName]
std::vector < std::pair <TString, TString> > getROOTobjects(TFile* f, TString rootDir, TString objectType);
std::vector < std::pair <TString, TString> > getROOTDirs(TFile* f, TString rootDir);


Bool_t map1DKeyExists(std::map <TString,TH1*> hist1D, TString histogram){
  /*  std::map <TString,TH1*>::const_iterator it = hist1D.find(histogram);
      return it! = hist1D.end();
  */
  return !(hist1D.find(histogram) == hist1D.end());
}

// Directory functions
Bool_t dirExists(TString dir);
void makeFullDir(TString baseDir, TString dir);
void makeDirectory(TString& directoryName);












//===================================================================
// Main function
//===================================================================


int getCalibration(){


  // Enable batch mode
  gROOT->SetBatch(kTRUE);


  unsigned int nParams = fitOffPT_vs_L1PT->GetNpar();

  // Create the directory in which to store the plots
  makeDirectory(plotDirectory);

  // list of histogram names in ROOT file with their corresponding directory structure
  std::vector < std::pair <TString, TString> > fileHistPair;
 
  // storage containers for histograms of various dimensions and their storage directory
  std::map <TString,TH1*> histogram1D;
  std::map <TString,TH2*> histogram2D;
  std::map <TString,TProfile*> histogramProf;
  std::map <TString,TString>   histogramDirectory;

  // temporary name of filepath and histogram
  TString filepath, filepathRaw, histogramName, histogramNameRaw;


  //================================================================
  // Automatically extract the histogram names from ROOT file and 
  // store them in maps
  //================================================================


  // open filestream
  TFile* f = new TFile(filename, "OPEN");

  // Check file exists
  if ( f->IsZombie() ) {
    std::cout << "\nERROR: File '" << filename << "' does not exist\n\n";
    exit(-1);
  }
  



  // Get a list of the main directories in the ROOT file
  std::vector < std::pair <TString, TString> > directories = getROOTobjects(f, ROOTdir, "TDirectoryFile");


  // Loop through all the directories
  for (uint iDir = 0; iDir < directories.size(); ++iDir){

    // Get current main directory
    ROOTdir = directories[iDir].second;
    std::cout << "Analysing the directory: " << ROOTdir << "\n\n\n";


    // *****************************************************************************************************************************
    // *                                                            TH2                                                            *
    // *****************************************************************************************************************************

#ifdef ANALYSE_TH2


    //=========================
    // 2D Histograms
    //=========================
    SUBPRINT("2D Histograms")
      // get list of histograms in file
      fileHistPair = getROOTobjects(f, ROOTdir, "TH2");


    // Subdirectories to run calibrations over
    std::vector<TString> subDirs;
//     subDirs.push_back( "/Calibration_PrePUS_ak5PUS_NVTXLt15/" );
//     subDirs.push_back( "/Calibration_PrePUS_ak5PUS_NVTXLt25/" );
//     subDirs.push_back( "/Calibration_PrePUS_ak5PUS_NVTXLt50/" );

    // ak5PUSRaw vs ak5PUS
    //    subDirs.push_back( "/Calibration_ak5PUSRaw_ak5PUS/" );

    subDirs.push_back( "/Calibration_PrePUS_ak5PUS_3Jets/" );

    //    subDirs.push_back( "/Calibration_PrePUS_ak5PUS_EtaLt3/" );
    //  subDirs.push_back( "/Calibration_PUS_ak5PUS/" );
    //  subDirs.push_back( "/Calibration_LPUS_ak5PUS/" );

    //	   subDirs.push_back( "/Calibration_UncalibCurr_ak5PUS/" );     


    std::map <TString, TString> typeLabel;
    typeLabel[ "/Calibration_ak5PUSRaw_ak5PUS/" ] = "ak5PUSRaw";

    typeLabel[ "/Calibration_PrePUS_ak5PUS_NVTXLt15/" ] = "PrePUS_NVTXLt15";
    typeLabel[ "/Calibration_PrePUS_ak5PUS_NVTXLt25/" ] = "PrePUS_NVTXLt25";
    typeLabel[ "/Calibration_PrePUS_ak5PUS_NVTXLt50/" ] = "PrePUS_NVTXLt50";

    typeLabel["/Calibration_PrePUS_ak5PUS/"]        = "PrePUS";
    typeLabel["/Calibration_PrePUS_ak5PUS_EtaLt3/"] = "PrePUSLt3";
    typeLabel["/Calibration_PUS_ak5PUS/"]           = "PUS";
    typeLabel["/Calibration_LPUS_ak5PUS/"]          = "LPUS";

    typeLabel["/Calibration_UncalibCurr_ak5PUS/"]   = "UncalibGCT";


    typeLabel[ "/Calibration_PrePUS_ak5PUS_3Jets/" ] = "PrePUS_3Jets";



    // Check labels are defined for each sample analysed
    for (uint iSub = 0; iSub < subDirs.size(); ++iSub ){

      // Check a label exists for the current subdir
      TString subDir = subDirs[ iSub ];
      std::map<TString, TString>::iterator labelIt = typeLabel.find( subDir );
      if ( labelIt == typeLabel.end() ){
	std::cout << "ERROR: No label was defined for subdir '" << subDir << "'\n\n"; 
	exit(1);
      }
      
    }


    // eta binning, key = lowerBin, value = binning range string. Used to obtain an eta ordered list of etabins
    std::map< double, TString > etaBins;

    for ( unsigned int iSub = 0; iSub < subDirs.size(); ++iSub ){

      // Get the current subdirectory to process
      TString curSubDir = subDirs[iSub];


      for (unsigned int i = 0;i < fileHistPair.size(); i++){

	// Extract the directory and name of the histogram
	filepathRaw      = fileHistPair[i].first;
	histogramNameRaw = fileHistPair[i].second;

	// pdf compatible filenames
	histogramName    = histogramNameRaw;
	filepath         = filepathRaw;

	// Find the directory with the calibration plots, currently e.g. : /JetHist/Calibration/Calibration_PrePUS_akPUS/iEtaBinned/iEta_-28to-25
	if (filepath.Contains("Calibration_") != 0){
	  // Only load the EtaBinned distributions
	  if (filepath.Contains("EtaBinned") != 0){
	    //	  if (filepath.Contains("iEtaBinned") != 0){

	    // Restrict to current subdirectory
	    if ( filepath.Contains( curSubDir ) != 0 ){
	      // TEMPORARY TO SPEED UP DEVELOPMENT RESTRICT TO RELEVENT COLLECTIONS
	      //	if ( (filepath.Contains("_PrePUS_") != 0) || (filepath.Contains("_PUS_") != 0) || (filepath.Contains("_LPUS_") != 0) ){
	      // Do not load NVTX binned distributions
	      if (filepath.Contains("NVTXBinned") == 0){
	  
	  
		// Restrict to calibration plots
		if ( ( histogramName.Contains("JetResponse_vs_OffPT") != 0 ) || ( histogramName.Contains("L1PT_vs_OffPT") != 0 ) ) {
			    
		  histogramName.ReplaceAll("(","{").ReplaceAll(")","}"); // Replace brackets to allow pdf printing
		  filepath.ReplaceAll("(","{").ReplaceAll(")","}");
	      
		  //	      std::cout << "NEW = " << filepath << "\t" << histogramName << "\n";
	      
		  // store histogram and its storage directory
		  histogram2D[ histogramName ]        = (TH2*)f->Get(filepathRaw + "/" + histogramNameRaw)->Clone();
		  histogramDirectory[ histogramName ] = filepath;
	      
		  // Ensure a directory is created to store the histogram
		  makeFullDir( plotDirectory, filepath );
	      
		  // change name to avoid memory conflicts
		  histogram2D[ histogramName ] ->SetName(histogramName);



		  // Get Eta label
		  TString EtaStr = histogramName;{
		  while ( EtaStr.Contains("Eta") )
		    EtaStr = EtaStr.Remove( 0, EtaStr.Index("Eta") + 4);
		  }
		  
		  TString etaLowStr  = EtaStr;
		  etaLowStr = etaLowStr.Remove( etaLowStr.Index("_to_") );
// 		  TString etaHighStr = EtaStr;
// 		  etaHighStr = etaHighStr.Remove( 0, etaHighStr.Index("_to_") + 4);  

		  double etaLow  = etaLowStr.Atof();
		  //		  double etaHigh = etaHighStr.Atof();
		  etaBins[ etaLow ] = EtaStr;


		}
	      }
	    }
	  }
	}
      }


      // Store an ordered vector of bins
      std::vector<TString> etaBinsOrdered;
      map<double, TString>::const_iterator itrBin;
      for (itrBin = etaBins.begin(); itrBin != etaBins.end(); ++itrBin){

	//	double etaLow  = itrBin->first;    // Extract bin lower eta
	TString etaBin = itrBin->second; 

	etaBinsOrdered.push_back( etaBin );

      }


   
      // ************************************************************************************
      // *                               Plotting histograms                                *
      // ************************************************************************************
      PRINT("Plotting 2D Histograms")
	// Enable the fit box
	gStyle->SetOptFit(111111);

      //=========================
      // 2D Histograms
      //=========================
      SUBPRINT("2D Histograms")
	// iterate through all histograms in the map
	map<TString, TH2*>::const_iterator itr2D;
  








      // **************************************************
      // *             Calibration parameters             *
      // **************************************************
  
      // pT width of bins
      int ptBinning = 20;
  
      // Range over which to make profiles
      double ptMin = 40; 
      //      double ptMin = 26; 

      //double ptMax = 120;
      double ptMax = 480;
      //      double ptMax = 140;
  
      // **************************************************
 
  
      // Determine the number of profile plots to produce (should probably round up?)
      int totalProfs = int(( ptMax - ptMin )/ptBinning);
  
  
      //DEBUGGING    
      //DEBUGGING
      //DEBUGGING
      //      totalProfs = 3;
      //DEBUGGING
      //DEBUGGING
      //DEBUGGING
  

      // Points stored for the calibration fit, indices = [iEta][ptBin]
      std::map<TString, std::vector< std::vector< std::pair<double,double> > > > collectioniEtaPTBinnedPoint;
      std::map<TString, std::vector< std::vector< std::pair<double,double> > > > collectioniEtaPTBinnedPointError;
  

  
//       //   // Points stored for the calibration fit, indices = [iEta][ptBin]
//       //   std::vector< std::vector< std::pair<double,double> > > iEtaPTBinnedPoint;
//       //   // Resize to fit the required data
//       //   iEtaPTBinnedPoint.resize( 56, std::vector< std::pair<double,double> >( totalProfs ) );

//       // Store whether to avoid processing the specified iEta
//       std::vector<bool> iEtaVeto(56, false);
  
  
      for(itr2D = histogram2D.begin(); itr2D != histogram2D.end(); ++itr2D){
    
	TString histName = itr2D->first;    // Extract histogram name
	TH2*    hist2D   = itr2D->second;   // Extract histogram
	filepath = histogramDirectory[ histName ] + "/"; // Extract directory
	TCanvas* canv = new TCanvas();




	// Get the current eta bin
	int etaIndex;
	TString EtaStr;
	for ( uint iEtaBin = 0; iEtaBin < etaBinsOrdered.size(); ++iEtaBin ){
	  
	  TString tempEtaBin = etaBinsOrdered[ iEtaBin ];
	  if ( histName.Contains( tempEtaBin ) != 0 ){
	    EtaStr = tempEtaBin;
	    etaIndex = iEtaBin;

	    std::cout << EtaStr << "\t" << histName << "\n";
	    break;
	  }

	}
	TString EtaStrLabel = EtaStr;
	EtaStrLabel.ReplaceAll("_to_", ", ");
	EtaStrLabel = "[" + EtaStrLabel + "]";



	hist2D->Draw("COLZ");                   
	canv->SaveAs(plotDirectory + filepath + histName + ".png");   // write histogram to file
	canv->SaveAs(plotDirectory + filepath + histName + ".pdf");   // write histogram to file

	// Store the name of the collection being calibrated, used to store the calibration data
	TString histBaseName = histName;

	if (histBaseName.Contains("JetResponse_vs_OffPT") != 0){
	  histBaseName = histBaseName.Remove( histBaseName.Index( "_JetResponse_vs_OffPT" ) );
	  //      std::cout << "RES " << histBaseName << "\n";
	}
	else if(histBaseName.Contains("L1PT_vs_OffPT") != 0 ){
	  histBaseName = histBaseName.Remove( histBaseName.Index( "_L1PT_vs_OffPT" ) );
	  //      std::cout << "PT " << histBaseName << "\n";
	}
	else{
	  std::cout << "ERROR: Could not extract the collection type from '" << histName << "'\n";
	  exit(1);
	}

	// Check whether storage exists for the current collection
	map<TString, std::vector< std::vector< std::pair<double,double> > > >::iterator it = collectioniEtaPTBinnedPoint.find( histBaseName );
	if ( it == collectioniEtaPTBinnedPoint.end() ) { 
	  // Doesn't exist, make new storage
      
	  // Points stored for the calibration fit, indices = [iEta][ptBin] 
	  std::vector< std::vector< std::pair<double,double> > > iEtaPTBinnedPoint;
	  std::vector< std::vector< std::pair<double,double> > > iEtaPTBinnedPointError;

	  int size = etaBinsOrdered.size();

	  // Resize to fit the required data 
	  iEtaPTBinnedPoint.resize( size, std::vector< std::pair<double,double> >( totalProfs ) );
	  collectioniEtaPTBinnedPoint[ histBaseName ] = iEtaPTBinnedPoint; 
	  iEtaPTBinnedPointError.resize( size, std::vector< std::pair<double,double> >( totalProfs ) );
	  collectioniEtaPTBinnedPointError[ histBaseName ] = iEtaPTBinnedPointError; 

	  std::cout << "Made collection '" << histBaseName << "'\n";
	}


	//DEBUGGING
	//DEBUGGING
	//DEBUGGING
	//     if (iEta != -28){
	//       continue;
	//     }
	//DEBUGGING
	//DEBUGGING
	//DEBUGGING


	// ****************************************************************************************************
	// *                                          Jet calibration                                         *
	// ****************************************************************************************************

	// Read details of the TH2
	int nBinsX      = hist2D->GetNbinsX();
	double binWidth = hist2D->GetBinWidth(1);
	double xPtMin   = hist2D->GetXaxis()->GetXmin();
    
	// Get current pt range to profile
	double ptLow, ptHigh;
    
	for ( int iBin = 0; iBin < totalProfs; ++iBin){
    
	  ptLow  = ptMin + iBin*ptBinning;
	  ptHigh = ptMin + (iBin + 1)*ptBinning - 1;
      
	  // Indices of the bins to project over, add one to account for underflow bin
	  int binLow  = (ptLow  - xPtMin)/binWidth + 1;
	  int binHigh = (ptHigh - xPtMin)/binWidth + 1;
      

	  //	std::cout << ptLow << " = " << hist2D->GetBinCenter( binLow ) << "\t" << ptHigh << " = " << hist2D->GetBinCenter( binHigh ) << "\n";

	  TString ptLowStr  =  TString(Form("%d",Int_t(ptLow)));
	  TString ptHighStr =  TString(Form("%d",Int_t(ptHigh)));
	  TString ptRangeStr = ptLowStr + "_to_" + ptHighStr;

	  TString label = typeLabel[ curSubDir ];

	  TString newSaveName = "_prof_" + ptRangeStr;
	  TString newHistName = "#eta #in " + EtaStrLabel + ", " + "p_{T} #in [" + ptLowStr + ", " + ptHighStr + "]";


	  if ( histName.Contains("JetResponse_vs_OffPT") != 0 ){
	    newHistName = label + " Response - " + newHistName;
	  }
	  else if ( histName.Contains("L1PT_vs_OffPT") != 0 ) {
	    newHistName = label + " L1 P_{T} - " + newHistName;
	  }


	  TH1D *yProject = hist2D->ProjectionY( newHistName, binLow, binHigh, "eo");




	  // ********************************************************************************
	  // *                           Jet response calibration                           *
	  // ********************************************************************************

    bool drawStats = false;
	  if ( histName.Contains("JetResponse_vs_OffPT") != 0 ){
     
	    //	    yProject->Rebin(5);
      if(!drawStats) yProject->SetStats(0);
      yProject->SetTitle(ptRangeStr+"GeV");
	    yProject->Draw();

	    TPaveStats *fitStat;

	    // ****************************************
	    // *             Fit gaussian             *
	    // ****************************************

	    double response    = -1;
	    double responseErr = -1;
      
      // Ensure the histogram contains entries
      if ( yProject->GetEntries() != 0){

        const double nSigma = 1.5;
        const double fitMin = -999;
        const int    nIter  = 3;

        fit_gaussian( yProject, nSigma, fitMin, nIter );



        TF1* gausResFit = (TF1*) yProject->GetListOfFunctions()->Last();

        //statsbox wizardry
        gPad->Update();
        if(drawStats){
          fitStat = (TPaveStats*)yProject->FindObject("stats");
          fitStat->SetTextColor(kBlue);
          fitStat->SetX1NDC(.57); fitStat->SetX2NDC(.88); fitStat->SetY1NDC(.14); fitStat->SetY2NDC(.53);
          fitStat->Draw();
        }

        if ( gausResFit) {
          //		yProject->Fit(gausResFit,"QR");
          gausResFit->SetLineColor(kRed);
          gausResFit->SetLineWidth(1);
          gausResFit->Draw("SAME");
        }
        canv->SaveAs(plotDirectory + filepath + histName + newSaveName + "FitGaus.png");   // write histogram to file
        canv->SaveAs(plotDirectory + filepath + histName + newSaveName + "FitGaus.pdf");   // write histogram to file


        if ( gausResFit) {

          // Get gaussian fit mean 
          response    = gausResFit->GetParameter( 1 );
          //responseErr = gausResFit->GetParameter( 2 );
          responseErr = gausResFit->GetParError( 1 );
        }	      
        std::cout << "Final fit: Response = " << response << " +/- " << responseErr << "\n";

      }




      std::cout << "Eta = " << EtaStr << "\tResponse = " << response << "\t+/- " << responseErr << "\n";	

      double reciprocalResponse = 0;
      double reciprocalResponseErr = 0;
      if (response != 0){
        reciprocalResponse    = 1/response;
        reciprocalResponseErr = responseErr/(response*response);
      }

      // Store the inverted mean for plotting
      collectioniEtaPTBinnedPoint[ histBaseName ][etaIndex][iBin].second      = reciprocalResponse;
      // Error in reciprocal is given by: responseErr / Response^2
      collectioniEtaPTBinnedPointError[ histBaseName ][etaIndex][iBin].second = reciprocalResponseErr;

    }



    // ********************************************************************************
    // *                              Jet PT calibration                              *
    // ********************************************************************************

    if ( (histName.Contains("L1PT_vs_OffPT") != 0) ){

      //	    yProject->Rebin(4);
      yProject->Draw();

      TPaveStats *fitStat;


      // ****************************************
      // *             Fit gaussian             *
      // ****************************************

      double l1Pt    = -1;
      double l1PtErr = -1;


      // Ensure the histogram contains entries
      if ( yProject->GetEntries() != 0){

        const double nSigma = 1.5;
        const double fitMin = -999;
        const int    nIter  = 3;


        fit_gaussian( yProject, nSigma, fitMin, nIter );


        TF1* gausPtFit = (TF1*) yProject->GetListOfFunctions()->Last();


        //statsbox wizardry
        gPad->Update();
        fitStat = (TPaveStats*)yProject->FindObject("stats");
        fitStat->SetTextColor(kBlue);
        fitStat->SetX1NDC(.57); fitStat->SetX2NDC(.88); fitStat->SetY1NDC(.14); fitStat->SetY2NDC(.53);
        fitStat->Draw();

        if ( gausPtFit) { 
          //		yProject->Fit(gausPtFit,"QR"); 
          gausPtFit->SetLineColor(kRed);
          gausPtFit->SetLineWidth(1);
          gausPtFit->Draw("SAME");
        }
        canv->SaveAs(plotDirectory + filepath + histName + newSaveName + "FitGaus.png");   // write histogram to file
        canv->SaveAs(plotDirectory + filepath + histName + newSaveName + "FitGaus.pdf");   // write histogram to file


        if (gausPtFit){
          // Get gaussian fit mean
          l1Pt    = gausPtFit->GetParameter( 1 );
          //l1PtErr = gausPtFit->GetParameter( 2 );
          l1PtErr = gausPtFit->GetParError( 1 );
        }

        std::cout << "Final fit: L1Pt = " << l1Pt << " +/- " << l1PtErr << "\n";

      }

      std::cout << "Eta = " << EtaStr << "\tl1Pt = " << l1Pt << "+/-\t" << l1PtErr << "\n";	
      // Store the mean for plotting
      collectioniEtaPTBinnedPoint[ histBaseName ][etaIndex][iBin].first      = l1Pt;
      collectioniEtaPTBinnedPointError[ histBaseName ][etaIndex][iBin].first = l1PtErr;

    }




  }


      }

      //#ifdef DEBUG_OFF  

      // ******************************************************************
      // *                    Plot calibration TGraphs                    *
      // ******************************************************************



      // E
      //      gStyle->SetOptFit(11111);




      std::map<TString, std::vector< std::vector< std::pair<double,double> > > >::const_iterator calibIt;

      for(calibIt = collectioniEtaPTBinnedPoint.begin(); calibIt != collectioniEtaPTBinnedPoint.end(); ++calibIt){


        TString collName = calibIt->first;    // Extract collection
        std::cout << "\n\nProcessing collection: " << collName << "------------------------------------------------------------\n\n\n\n";


        // ROOT file to store calibration TGraphs
        TFile* graphFile = new TFile( plotDirectory + collName + ".root", "RECREATE");

        // Extract the directory to store plots
        TString calibFilePath = filepath;
        calibFilePath = calibFilePath.Remove( calibFilePath.Index( "EtaBinned" ) );


        // Create a CMSSW compatable LUT
        TString outputCMSSW = "\t" + collName + "LUT = cms.vdouble(\n";

        // Iterate over iEta bins
        //	for (int iEtaIndex = 0; iEtaIndex < 56; ++iEtaIndex){
        for (uint etaIndex = 0; etaIndex < etaBinsOrdered.size(); ++etaIndex){

          TString EtaStr = etaBinsOrdered[ etaIndex ];
          TString EtaStrLabel = EtaStr;
          EtaStrLabel.ReplaceAll("_to_", ", ");
          EtaStrLabel = "[" + EtaStrLabel + "]";

          // Make a TGraph for each iEta bin
          double xArr[999],    yArr[999];
          double xErrArr[999], yErrArr[999];
          // Residuals
          double xResArr[999], yResArr[999];

          for (int iProf = 0; iProf < totalProfs; ++iProf){

            double x    = collectioniEtaPTBinnedPoint[collName][etaIndex][iProf].first;
            double y    = collectioniEtaPTBinnedPoint[collName][etaIndex][iProf].second;
            double xErr = collectioniEtaPTBinnedPointError[collName][etaIndex][iProf].first;
            double yErr = collectioniEtaPTBinnedPointError[collName][etaIndex][iProf].second;

            // Skip failed fits
            if ( ( x == -1) || ( y == -1) ){ continue; }

            std::cout << "[ " << x << ", " << y << " )\n";

            // Store data in arrays for TGraph input
            xArr[ iProf ] = x;
            yArr[ iProf ] = y;
            xErrArr[ iProf ] = xErr;
            yErrArr[ iProf ] = yErr;

          }

          // Make the Eta binned TGraph
          TGraphErrors*  calibGraph = new TGraphErrors( totalProfs, xArr, yArr, xErrArr, yErrArr );
          calibGraph->SetTitle("#eta #in " + EtaStrLabel + ";<L1 p_{T}> (GeV);<L1 p_{T}/RECO p_{T}>^{-1}");
          calibGraph->SetMarkerStyle(8);
          calibGraph->SetMarkerSize(0.5);

          TCanvas* canv = new TCanvas();
          canv->SetGridx();
          canv->SetGridy();

          calibGraph->Draw("AP*");
          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_iEta_" + EtaStr + ".png");   // write histogram to file    
          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_iEta_" + EtaStr + ".pdf");   // write histogram to file    

          TString EtaStrROOT = EtaStr;
          EtaStrROOT.ReplaceAll("-","minus");

          // Save to the ROOT file
          calibGraph->Write( "Eta_" + EtaStrROOT );



          // ********************************************************************************
          //                                     FITTING                                    *
          // ********************************************************************************

#ifdef FITTING

          double calibFitMin = 20;
          double calibFitMax = 230;

          // Fit the calibration curve over a restricted range




          // ********************************************************************************
          // FROM: jet_l3_correction_x.cc

          // response 
          TF1* fitrsp;
          fitrsp->SetLineColor(kRed);
          fitrsp->SetLineWidth(1);

          fitrsp = new TF1("fitrsp","[0]-[1]/(pow(log10(x),2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))",
              1.0,calibGraph->GetX()[calibGraph->GetN()-1]);

          fitrsp->SetParameter(0,0.96);
          fitrsp->SetParameter(1,0.033);
          fitrsp->SetParameter(2,-0.7);
          fitrsp->SetParameter(3,0.02);
          fitrsp->SetParameter(4,1.02);
          fitrsp->SetParameter(5,2.7);
          fitrsp->SetParameter(6,0.016);

          calibGraph->Fit(fitrsp,"QR");

          //statsbox wizardry
          TPaveStats *fitStat;
          gPad->Update();
          fitStat = (TPaveStats*)calibGraph->FindObject("stats");
          fitStat->SetTextColor(kRed);
          fitStat->Draw();

          // ********************************************************************************

          // ********************************************************************************
          // FROM: jet_l3_correction_x.cc

          // 	  // L3Absolute correction 
          // 	  string fitcor_as_str;

          // 	  fitcor_as_str = "[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";


          // 	  TF1* fitcor = new TF1("fitcor",fitcor_as_str.c_str(),
          // 				2.0,gcor->GetX()[gcor->GetN()-1]);

          // 	  (alg.find("pf")!=string::npos) {
          // 	    fitcor->SetParameter(0,1.04);
          // 	    fitcor->SetParameter(1,.033);
          // 	    fitcor->SetParameter(2,-0.7);
          // 	    fitcor->SetParameter(3,0.02);
          // 	    fitcor->SetParameter(4,1.02);
          // 	    fitcor->SetParameter(5,2.7);
          // 	    //  fitrsp->SetParameter(6,0.016);  // ? BUG ?

          // 	    gcor->Fit(fitcor,"QR");

          // ********************************************************************************


          // 	  TF1 *calibrationFit = new TF1("calibrationFit", calibFit, calibFitMin, calibFitMax, 6); 
          // 	  calibrationFit->SetLineColor(kRed);
          // 	  calibrationFit->SetLineWidth(1);

          // 	  // Set initial parameter values
          // 	  calibrationFit->SetParameter(0, 0.78);                                                   
          // 	  calibrationFit->SetParameter(1, 4.33);                                                   
          // 	  calibrationFit->SetParameter(2, 2.67);                                                   
          // 	  calibrationFit->SetParameter(3, 0.57);                                                   
          // 	  calibrationFit->SetParameter(4, 0.88);                                                   
          // 	  calibrationFit->SetParameter(5, 0.41);                                                   

          // 	  calibGraph->Fit( calibrationFit, "MR");

          // 	  //statsbox wizardry
          // 	  TPaveStats *fitStat;
          // 	  gPad->Update();
          // 	  fitStat = (TPaveStats*)calibGraph->FindObject("stats");
          // 	  fitStat->SetTextColor(kRed);
          // 	  //	fitStat->SetX1NDC(.27); fitStat->SetX2NDC(.58); fitStat->SetY1NDC(.54); fitStat->SetY2NDC(.93);
          // 	  fitStat->Draw();


          //     // Fit the calibration curve over the entire range
          //     TF1 *calibrationFitFull = new TF1("calibrationFit",calibFit,5,250,6); 
          //     calibrationFitFull->SetLineStyle(7);
          //     calibrationFitFull->SetLineWidth(1);
          //     // Set initial parameter values
          //       calibrationFitFull->SetParameter(0, 0.78);                                                   
          //       calibrationFitFull->SetParameter(1, 4.33);                                                   
          //       calibrationFitFull->SetParameter(2, 2.67);                                                   
          //       calibrationFitFull->SetParameter(3, 0.57);                                                   
          //       calibrationFitFull->SetParameter(4, 0.88);                                                   
          //       calibrationFitFull->SetParameter(5, 0.41);                                                   

          //     calibGraph->Fit( calibrationFitFull, "+MR");



          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_Eta_" + EtaStr + "_FIT.png");   // write histogram to file 
          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_Eta_" + EtaStr + "_FIT.pdf");   // write histogram to file 




          // Save to the ROOT file
          calibGraph->Write( "Eta_" + EtaStrROOT );


          // Caculate the fit residuals
          for (int iProf = 0; iProf < totalProfs; ++iProf){


            double x = xArr[ iProf ];
            double y = yArr[ iProf ];

            double fitY = calibrationFit->Eval(x);

            xResArr[iProf] = x;
            yResArr[iProf] = fitY - y;

          }


          // Plot the residual of the fit
          TGraph* calibRes = new TGraph( totalProfs, xResArr, yResArr );
          calibRes->SetTitle("#eta #in " + EtaStrLabel + ";<L1 p_{T}> (GeV);Residual <L1 p_{T}/RECO p_{T}>^{-1}");
          calibRes->SetMarkerStyle(8);
          calibRes->SetMarkerSize(0.5);
          calibRes->Draw("AP*");

          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_Eta_" + EtaStr + "_FIT_Res.png");   // write histogram to file 
          canv->SaveAs(plotDirectory + calibFilePath + "calibGraph_Eta_" + EtaStr + "_FIT_Res.pdf");   // write histogram to file 
          std::cout << "\n\nPLINT: " << plotDirectory + calibFilePath + "calibGraph_Eta_" + EtaStr + "_FIT_Res.png" << "\n\n";

          calibRes->Write( "Eta_" + EtaStrROOT + "_Res" );



          // Extract calibration parameters
          TString p0 = Form("%f", calibrationFit->GetParameter( 0 ));
          TString p1 = Form("%f", calibrationFit->GetParameter( 1 ));
          TString p2 = Form("%f", calibrationFit->GetParameter( 2 ));
          TString p3 = Form("%f", calibrationFit->GetParameter( 3 ));
          TString p4 = Form("%f", calibrationFit->GetParameter( 4 ));
          TString p5 = Form("%f", calibrationFit->GetParameter( 5 ));

          // Output fit parameters
          std::cout << "\tp0 = " << p0 << "\n"
            << "\tp1 = " << p1 << "\n"
            << "\tp2 = " << p2 << "\n"
            << "\tp3 = " << p3 << "\n"
            << "\tp4 = " << p4 << "\n"
            << "\tp5 = " << p5 << "\n\n\n";

          outputCMSSW += "\t\t" + p0 + ",\t" + p1 + ",\t" + p2 + ",\t" + p3 + ",\t" + p4 + ",\t" + p5 + ",\n";

#endif

        }

        outputCMSSW += "\n\t),";

        ofstream calibLUT;
        calibLUT.open( plotDirectory + collName + "LUT.txt" );
        calibLUT << outputCMSSW;
        calibLUT.close();


        std::cout << "Made LUT: " + plotDirectory + collName + "LUT.txt" + "\n";
        std::cout << outputCMSSW << "\n\n";


        std::cout << "\n\nMade ROOT file: " + plotDirectory + collName + ".root" + "\n";
        // Close the ROOT file
        graphFile->Close();

        }

        //#endif

        // Delete the previous contents of the container
        histogram2D.clear();

      }

#endif








      // Delete the previous contents of the container
      histogramDirectory.clear();


      } // Directory loop







      // **************************************************************************************************
      // *                                          Provenance                                            *
      // **************************************************************************************************


      // To test different text positions run root in interactive mode with the code:
      // TCanvas canv;
      // canv.SetGrid();
      // canv.DrawFrame(0,0,1,1);
      //
      // Print statistics of the program running
      TCanvas statCanv;
      TLatex latex;
      Float_t curY = 0.95; // current Y coordinate of text
      latex.SetTextSize(0.04);
      latex.DrawLatex(0.05, curY, "Input file:");
      latex.SetTextSize(0.025);
      curY -= 0.02;;
      latex.DrawLatex(0.05, curY, filename);


      // save and close .eps
      statCanv.SaveAs( plotDirectory + "Provenance.png" );
      statCanv.SaveAs( plotDirectory + "Provenance.pdf" );



      // Make webpage
      chdir( plotDirectory.Data() );
      system( "/home/hep/mb1512/.scripts/makeHTML.py" );


      std::cout << "\n\n\nOutput histograms to webpage: \n\n\t"
        << plotDirectory.ReplaceAll("/home/hep/", "http://www.hep.ph.ic.ac.uk/~").ReplaceAll("/public_html","") 
        << "\n\n\n";//mb1512/public_html/plots/";

      f->Close();
      return 0;

  }

  // --------------------------------------------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------------------------------------------








  // // Stores all histograms in memory which are of a given TObject type

  // void loadHistograms(TFile* f, TString rootDir, TString histType, std::vector < pair <TString, TString> > fileHistPair,
  // 		    std::map <TString,TH1*> histMap){


  //   // get list of histograms in file of specified type
  //   fileHistPair = getROOTobjects(f, rootDir, histType);

  //   for (unsigned int i = 0;i < fileHistPair.size(); i++){

  //     // Extract the directory and name of the histogram
  //     filepath         = fileHistPair[i].first;
  //     histogramNameRaw = fileHistPair[i].second;

  //     // create a .pdf compatible filename
  //     histogramName    = histogramNameRaw;
  //     histogramName.ReplaceAll("(","{").ReplaceAll(")","}"); // Replace brackets to allow pdf printing
  //     filepath.ReplaceAll("(","{").ReplaceAll(")","}");

  //     std::cout << "NEW = " << filepath << "\t" << histogramName << "\n";


  //     // store histogram and its storage directory
  //     histogramMap[ histogramName ]       = (TH1*)f->Get(filepath + "/" + histogramNameRaw)->Clone();
  //     histogramDirectory[ histogramName ] = filepath;

  //     // Ensure a directory is created to store the histogram
  //     makeFullDir( plotDirectory, filepath );

  //     // change name to avoid memory conflicts
  //     histogramMap[ histogramName ] ->SetName(histogramName);

  //   }

  // }














  //##############################################
  //## TFile functions
  //##############################################

  //Returns a boolean value of whether the object in the specified directory exists
  //Checks first whether directory exists to avoid Segfaults
  Bool_t objectExists(TFile* file,  TString directory, TString objectName){
    //Key to determine whether specified objects exist
    TKey *key;

    if (directory != ""){

      //Check if the directory has more than one level
      if (directory.Contains("/")){
        std::vector <TString> dirLevel = splitString(directory,"/");
        TString curDir = "";

        //Loop through all the levels of the directory
        for (unsigned int i = 0;i < dirLevel.size();i++){
          if (curDir == "")
            //Check first level exists
            key = file->FindKey(dirLevel[i]);
          else{
            //Check higher levels exist
            key = file->GetDirectory(curDir)->FindKey(dirLevel[i]);
          }

          //Check if the directory was not found
          if (key == 0)	
            return false;

          curDir += dirLevel[i] + "/";
        }
      }
      else{
        //Check directory exists for a single level directory
        key = file->FindKey(directory);
      }
    }


    //Directory exists so now attempt to find the object
    key = file->GetDirectory(directory)->FindKey(objectName);

    //Return false if the object was not found
    return (key != 0);
  }
  //##############################################


  std::vector <TString> splitString(TString inputStr, TString splitChar){

    //Vector to store the split substrings
    std::vector <TString> splitArr;

    //Initialise values so that it is detectable if the split character is not present
    Int_t index(-1), tempIndex(-1);

    //Find the first index of the split character, returns -1 if char not found
    tempIndex = inputStr.Index(splitChar,0);	

    TString tempStr, inputCopy;

    //Keep looping while the character is found
    while(tempIndex != -1){

      //KEEP THE SUBSRING LEFT OF THE STRIP CHARACTER
      inputCopy = inputStr; //Make a fresh copy of the input string
      inputCopy.Remove(tempIndex,inputCopy.Length());
      //if (index != 0)
      //  index +=1;
      inputCopy.Remove(0,index+1);

      //Store if the substring isn't empty
      if (inputCopy != "")
        splitArr.push_back(inputCopy);

      index = tempIndex;
      //Search for next instance of the character from index tempIndex+1 onwards
      tempIndex = inputStr.Index(splitChar,tempIndex+1);

    }

    //Check that the split character was found
    if (index == -1){
      std::cout << "\nError the character: \'" << splitChar 
        << "\' was not found in the string: \"" << inputStr << "\"\n\n";
      exit(1);
    }

    //Recover the last substring which is right of the final character

    inputCopy = inputStr;
    tempStr = inputCopy.Remove(0,index+1);

    if (tempStr != "")
      splitArr.push_back(tempStr);


    return splitArr;
  }



  // Returns a list of the specified objectType in a given rootDir
  // "TDirectoryFile"
  // [ntuplePath, objectName]
  std::vector < std::pair <TString, TString> > getROOTobjects(TFile* f, TString rootDir, TString objectType){

    std::vector < std::pair <TString, TString> > objectList;

    TIter nextkey(f->GetDirectory(rootDir)->GetListOfKeys());
    TKey *key;
    TString objName, objClassName;
    TString ntuplePath;

    while ((key=(TKey*)nextkey())) {  //Loop while the current key is valid

      TObject *obj = key->ReadObj(); //Point to the object with the current key

      objName      = key->GetName();
      objClassName = obj->ClassName();

      //    if ( objClassName == objectType ){
      if ( objClassName.Contains(objectType) ){

        ntuplePath = rootDir;
        objectList.push_back( std::make_pair(rootDir,objName) );
        //      std::cout << objName << "\t" << objClassName << "\n";

      }
      else if ( objClassName == "TDirectoryFile"){


        //       std::cout << "----------------------------------------------------------------------------------\n\n";
        //       std::cout << objName << "\t" << objClassName << "\n\n";

        // Object is a directory so ammend the new ntuple path
        ntuplePath = rootDir + "/" + objName;
        // Look for the specified object in this new ntuple path
        std::vector < std::pair <TString, TString> > subDirObjectList = getROOTobjects(f, ntuplePath, objectType);

        if (subDirObjectList.size() == 0){
          //	std::cout << "No objects of type " << objectType << " found in directory: " << ntuplePath << "\n";
        }
        else{

          // Add the found objects to the stored list
          objectList.insert(objectList.end(), subDirObjectList.begin(), subDirObjectList.end());

        }


      }

    }



    return objectList;
    }



    std::vector < std::pair <TString, TString> > getROOTDirs(TFile* f, TString rootDir){

      std::vector < std::pair <TString, TString> > objectList;

      TIter nextkey(f->GetDirectory(rootDir)->GetListOfKeys());
      TKey *key;
      TString objName, objClassName;
      TString ntuplePath;

      TString objectType = "TDirectoryFile";


      while ((key=(TKey*)nextkey())) {  //Loop while the current key is valid

        TObject *obj = key->ReadObj(); //Point to the object with the current key

        objName      = key->GetName();
        objClassName = obj->ClassName();


        if ( objClassName.Contains(objectType) ){

          ntuplePath = rootDir;
          objectList.push_back( std::make_pair(rootDir,objName) );
          std::cout << objName << "\t" << objClassName << "\n";

          // Object is a directory so ammend the new ntuple path                                                                                                   
          ntuplePath = rootDir + "/" + objName;
          // Look for the specified object in this new ntuple path
          std::vector < std::pair <TString, TString> > subDirObjectList = getROOTDirs(f, ntuplePath);

          if ( !(subDirObjectList.empty()) ){

            // Add the found objects to the stored list 
            objectList.insert(objectList.end(), subDirObjectList.begin(), subDirObjectList.end());

          }

        }

      }



      return objectList;
    }





    //##############################################
    //## Directory functions
    //##############################################

    // Returns a boolean value of whether the directory exists
    Bool_t dirExists(TString dir){
      //Check if directory already exists
      struct stat sb; 

      if (stat(dir, &sb) == 0) // stat == 0 <=> Directory exists
        return true;
      else
        return false;
    }


    // Makes the entire file structure required to create the uppermost directory
    void makeFullDir(TString baseDir, TString dir){


      if ( !(dir.Contains("/")) ){
        dir += "/";
      }

      // Generate list of all the subdirectories
      std::vector <TString> subDirs = splitString(dir,"/");

      TString subDir = "";

      // Check each directory exists to the uppermost directory
      for (unsigned int iSplit = 0; iSplit < subDirs.size(); iSplit++){


        TString curDir = subDirs[iSplit];
        subDir += curDir + "/";

        TString filePath = baseDir + subDir;

        //      std::cout << "LAME    " << filePath << "\n";


        if ( !dirExists(filePath)){
          //Make the directory with the corresponding RW permissions
          mkdir( filePath.Data(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
          //      std::cout << "Made directory: " << filePath << "\n";
        }

      }


    }

    // Takes a directory name and makes a unique directory with a revision number.
    // Modifies the input directory so that it matches this new directory name.
    void makeDirectory(TString& directoryName){


      // Get current date to use for a folder name
      time_t rawTime;
      struct tm* timeInfo;
      char curDate[80];

      time(&rawTime);
      timeInfo = localtime(&rawTime);
      strftime(curDate,80,"%d_%m_%y",timeInfo);  // Extract date in format : dd_mm_yyyy

      //Store the directory name with date
      directoryName.Append(curDate);

      //Check if directory already exists, if so add a 'revision' prefix
      if (dirExists(directoryName)){
        Int_t num = 1;
        TString tempStr = "";

        do{

          tempStr = directoryName + "_rev";
          tempStr += num;

          num++;
        }while(dirExists(tempStr)); //Loop while the pathname already exists

        //Store the directory filename
        directoryName = tempStr;
      }

      directoryName += "/";

      //Make the directory with the corresponding RW permissions
      mkdir(directoryName.Data(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    }




    //______________________________________________________________________________ 
    //
    // JETMET additions
    //______________________________________________________________________________                                                                              
    void fit_gaussian( TH1D*& hrsp, const double nsigma, const double fitMin, const int niter){ 

      if (0==hrsp) {
        cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
      }

      string histname = hrsp->GetName();
      double mean     = hrsp->GetMean();
      double rms      = hrsp->GetRMS();
      double ptRefMax(1.0),rspMax(0.0);

      double norm  = hrsp->GetMaximumStored();
      double peak  = mean;
      double sigma = rms;

      std::cout << "\n\nInitial fit: Norm = " << norm << "\tPeak = " << peak << "\tSigma = " << sigma << "\n"; 

      //   int pos1     = histname.find("RefPt");
      //   int pos2     = histname.find("to",pos1);
      //   string ss    = histname.substr(pos1+5,pos2);
      //   if (from_string(ptRefMax,ss,std::dec)) {
      //     if (histname.find("RelRsp")==0)
      //       rspMax = jtptmin/ptRefMax;
      //     if (histname.find("AbsRsp")==0)
      //       rspMax = jtptmin-ptRefMax;
      //   }


      // Fit range
      double xmin  = hrsp->GetXaxis()->GetXmin();
      double xmax  = hrsp->GetXaxis()->GetXmax();
      TF1* fitfnc(0); int fitstatus(-1);

      // Perform the specified number of fitting iterations
      for (int iiter=0;iiter<niter;iiter++) {

        vector<double> vv;
        vv.push_back( fitMin );
        vv.push_back(xmin);
        vv.push_back(peak-nsigma*sigma);
        // Modify fit range
        double fitrange_min = *std::max_element(vv.begin(),vv.end());
        double fitrange_max = std::min(xmax,peak+nsigma*sigma);

        std::cout << "Fit range: " << fitrange_min << " " << fitrange_max << std::endl;
        adjust_fitrange( hrsp, fitrange_min, fitrange_max );



        // Create new function fitting over specified range
        fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
        fitfnc->SetParNames("N","#mu","#sigma");
        fitfnc->SetParameter(0,norm);
        fitfnc->SetParameter(1,peak);
        fitfnc->SetParameter(2,sigma);
        fitstatus = hrsp->Fit(fitfnc,"RQ0");
        delete fitfnc;
        fitfnc = hrsp->GetFunction("fgaus");
        //fitfnc->ResetBit(TF1::kNotDraw); 
        if (fitfnc) {
          norm  = fitfnc->GetParameter(0);
          peak  = fitfnc->GetParameter(1);
          sigma = fitfnc->GetParameter(2);
          std::cout << "Fit iteration = " << iiter << "\tNorm = " << norm << "\tPeak = " << peak << "\tSigma = " << sigma << "\n"; 
        }
      }
      if(hrsp->GetFunction("fgaus")==0)
      {
        cout << "No function recorded in histogram " << hrsp->GetName() << endl;
      }
      if (0!=fitstatus){
        cout<<"fit_gaussian() to "<<hrsp->GetName()
          <<" failed. Fitstatus: "<<fitstatus
          <<" - FNC deleted."<<endl;
        hrsp->GetListOfFunctions()->Delete();
      }
    }

    //______________________________________________________________________________ 

    void adjust_fitrange(TH1* h,double& min,double& max)
    {
      int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
      int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
      while ((imax-imin)<8) {
        if (imin>1) {imin--; min = h->GetBinCenter(imin); }
        if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
      }
    }

