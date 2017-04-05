#include <iostream>
#include <fstream>

//#include "TFile.h"
//#include "TTree.h"

using namespace std;

int main (int argc, char** argv)
{
  string steeringFileName = argv[1];
  
  ifstream steeringFile;
  steeringFile.open(steeringFileName);
  
  string word;
  string filename;
  double pressure;
  double timeRes;
  
  while (steeringFile >> word)
  {
    filename = word;
    
    steeringFile >> word;
    pressure << atof(word);
    
    steeringFile >> word;
    timeRes << atof(word);
    
    cout << "Reading file: " << filename << " (pressure " << pressure << " atm, time resolution " << timeRes << " ns)" << endl;
  }
  
  steeringFile.close();
  
  return 0;
}
