#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1I.h"
#include "TGraph2D.h"

using namespace std;

int main (int argc, char** argv)
{
  string steeringFileName = argv[1];
  
  ifstream steeringFile;
  steeringFile.open(steeringFileName.c_str());
  
  string word;
  string filename;
  Double_t pressure;
  Double_t timeRes;
  
  TFile* analysisFile = new TFile("/storage/gp_ws_ddm/Analysis.root", "UPDATE");
  
  TGraph2D* deviationPlot = new TGraph2D(1);
  
  Int_t pointCounter = 0;
  
  while (steeringFile >> word)
  {
    pointCounter++;
    deviationPlot->Set(pointCounter);
    
    TH1I* deviationHist = new TH1I("deviation", "deviation", 100, 0, M_PI);
    
    // Get parameters from steering file
    
    steeringFile >> word;
    pressure = atof(word.c_str());
    
    steeringFile >> word;
    timeRes = atof(word.c_str());
    
    stringstream filename;
    filename << "/storage/gp_ws_ddm/simOutput/ddm_p" << pressure << "atm_tres" << timeRes << "ns.root";
    
    cout << "Reading file: " << filename.str() << " (pressure " << pressure << " atm, time resolution " << timeRes << " ns)" << endl;
    
    // Open data file and read in results tree
    TFile* dataFile = new TFile(filename.str().c_str(), "READ");
    TTree* dataTree = (TTree*)dataFile->Get("recoResultsCamera");
    TBranch* dataBranch = dataTree->GetBranch("recoResultsCamera_branch");
    Int_t nEvents = dataTree->GetEntries();
    
    // Get branch contents
    Double_t branchData[9];
    dataTree->SetBranchAddress("recoResultsCamera_branch", &branchData);
    
    // Loop over data
    for (Int_t i=0; i<nEvents; i++)
    {
      dataTree->GetEntry(i);
      
      deviationHist->Fill(branchData[6]);
    }
    
    // Calculate mean and median of deviation distribution
    Double_t deviationMean = deviationHist->GetMean();
    Double_t deviationMedian;
    Double_t medianQuantile = 0.5;
    
    deviationHist->GetQuantiles(1,&deviationMedian,&medianQuantile);
    
    cout << "Median deviation: " << deviationMedian << endl << endl;
    
    delete deviationHist;
    
    // Close data file
    dataFile->Close();
    
    // Save data point
    deviationPlot->SetPoint(pointCounter - 1, pressure, timeRes, deviationMedian);
  }
  
  deviationPlot->Write();
  delete deviationPlot;
  
  steeringFile.close();
  analysisFile->Close();
  
  return 0;
}
