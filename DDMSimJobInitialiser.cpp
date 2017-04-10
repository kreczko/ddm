#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

int main()
{
  ofstream shellScript;
  ofstream steeringFile;
  
  int numPointsP = 10; // number of points (pressure)
  int numPointsTres = 10; // number of points (time resolution)
  
  double lowP = 0.005; // lowest pressure (atm)
  double lowTres = 1.0; // lowest time resolution (ns)
  
  double stepP = 0.005; // step in pressure (atm)
  //double stepTres = 10.0; // step in time resolution (ns)
  double factorTres = 3.0; // multiplicative factor for time resolution (ns)
  
  shellScript.open("/storage/gp_ws_ddm/DDMSimJob.sh");
  steeringFile.open("/storage/gp_ws_ddm/simOutput/steering.txt", fstream::app);
  
  shellScript << "# Pressure (atm): " << lowP << "--" << lowP + ((numPointsP-1)*stepP) << endl;
  //shellScript << "# Time resolution (ns): " << lowTres << "--" << lowTres + ((numPointsTres-1)*stepTres) << endl;
  shellScript << "# Time resolution (ns): " << lowTres << "--" << lowTres * pow(factorTres, numPointsTres-1) << endl;
  shellScript << "cd ~/geant4/DDM-build" << endl;
  
  int seed = 2;
  
  for (int i = 0; i < numPointsP; i++)
  {
    double pressure = lowP + (i*stepP);
      
    for (int j = 0; j < numPointsTres; j++)
    {
      double timeRes = lowTres * pow(factorTres,j);
      
      ofstream paramFile;
      
      stringstream filename;
      filename << "/storage/gp_ws_ddm/parameterFiles/P" << pressure << "atm_Tres" << timeRes << "ns.txt";
      
      paramFile.open(filename.str().c_str());
      paramFile << "Pressure: " << pressure << endl;
      paramFile << "TimeRes: " << timeRes << endl;
      paramFile.close();
      
      shellScript << "./DDM -s 1 -storage 1 -m isoArgon.mac -p " << filename.str() << " -r " << seed << endl;
      steeringFile << filename.str() << " " << pressure << " " << timeRes << endl;
      
      seed++;
    }
  }
  
  shellScript << "echo === SIMULATION JOB COMPLETE ===";
  
  shellScript.close();
  steeringFile.close();
  
  return 0;
}
