#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

int main()
{
  ofstream shellScript;
  
  int numPointsP = 2; // number of points (pressure)
  int numPointsTres = 2; // number of points (time resolution)
  
  double lowP = 0.005; // lowest pressure (atm)
  double lowTres = 10.0; // lowest time resolution (ns)
  
  double stepP = 0.005; // step in pressure (atm)
  double stepTres = 10.0; // step in time resolution (ns)
  
  shellScript.open("/storage/gp_ws_ddm/DDMSimJob.txt");
  
  shellScript << "# Pressure (atm): " << lowP << "--" << lowP + (numPointsP*stepP) << endl;
  shellScript << "# Time resolution (ns): " << lowTres << "--" << lowTres + (numPointsTres*stepTres) << endl;
  shellScript << "cd ~/geant4/DDM-build" << endl;
  
  for (int i = 0; i < numPointsP; i++)
  {
    double pressure = lowP + (i*stepP);
      
    for (int j = 0; j < numPointsTres; j++)
    {
      double timeRes = lowTres + (j*stepTres);
      
      ofstream paramFile;
      
      stringstream filename;
      filename << "/storage/gp_ws_ddm/parameterFiles/P" << pressure << "atm_Tres" << timeRes << "ns.txt";
      
      paramFile.open(filename.str().c_str());
      paramFile << "Pressure: " << pressure << endl;
      paramFile << "TimeRes: " << timeRes << endl;
      paramFile.close();
      
      shellScript << "./DDM -s 1 -storage 1 -m isoArgon.mac -p " << filename << endl;
    }
  }
  
  shellScript << "echo === SIMULATION JOB COMPLETE ===";
  
  shellScript.close();
  
  return 0;
}
