#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"

using namespace std;

int main (int argc, char** argv)
{
  string steeringFileName = argv[1];
  
  ifstream steeringFile;
  steeringFile.open(steeringFileName.c_str());
  
  string word;
  string filename;
  double pressure;
  double timeRes;
  
  TFile* analysisFile = new TFile("/storage/gp_ws_ddm/Analysis.root", "UPDATE");
  
  while (steeringFile >> word)
  {
    // Get file name and parameters from steering file
    filename = word;
    
    steeringFile >> word;
    pressure = atof(word.c_str());
    
    steeringFile >> word;
    timeRes = atof(word.c_str());
    
    cout << "Reading file: " << filename << " (pressure " << pressure << " atm, time resolution " << timeRes << " ns)" << endl;
    
    // Open data file
    TFile* analysisFile = new TFile(filename.c_str(), "READ");
    
    // Close data file
    analysisFile->Close();
    
  }
  
  steeringFile.close();
  analysisFile->Close();
  
  return 0;
}
