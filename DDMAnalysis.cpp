#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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
    // Get parameters from steering file
    
    
    steeringFile >> word;
    pressure = atof(word.c_str());
    
    steeringFile >> word;
    timeRes = atof(word.c_str());
    
    stringstream filename;
    filename << "/storage/gp_ws_ddm/simOutput/ddm_p" << pressure << "atm_tres" << timeRes << "ns.root";
    
    cout << "Reading file: " << filename.str() << " (pressure " << pressure << " atm, time resolution " << timeRes << " ns)" << endl;
    
    // Open data file
    TFile* analysisFile = new TFile(filename.str().c_str(), "READ");
    
    // Close data file
    analysisFile->Close();
    
  }
  
  steeringFile.close();
  analysisFile->Close();
  
  return 0;
}
