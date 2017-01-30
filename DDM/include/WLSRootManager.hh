#ifndef ROOT_MANAGER_HH
#define ROOT_MANAGER_HH

#include "TFile.h"
#include "TTree.h"
#include "globals.hh"
#include <string>
#include "G4String.hh"
#include "TH1D.h"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <string>

using namespace std;

class WLSRootManager
{
	private:
		TFile* root_file;
		//Define _mng variables/arrays as needed.
		G4double PhotonInformation_mng[9];
		G4double BoxInformation_mng[4];
		G4double MacroPhotonData_mng[9];
		string ClockTime_mng;
		G4String GpsData_mng;

		G4double EventID_mng;
		G4double RunID_mng;

		G4int totalPhotons_mng;
		G4int photonDetectCounter_mng;
		G4int photonAbsorbCounter_mng;
		G4int photonEscapeCounter_mng;

		G4int photonBounceCounter_mng;
		G4int overallBounceCounter_mng;

		G4double overallEnergyCounter_mng;
		G4double overallLengthCounter_mng;
		G4double overallLifetimeCounter_mng;

		G4bool batchSpecialInputs_mng;
		G4double batchGeomX_mng;
		G4double batchGeomY_mng;
		G4double batchGeomZ_mng;

	public:
		static void CreateRootManager(G4String filename);
		static WLSRootManager* GetRootManager();
		static void DestroyRootManager();

		WLSRootManager(G4String filename);

		void InitialiseTree(G4String treename1, G4String treename2, G4String treename3);//Add more trees as arguments if desired.

		//void FillTree(); Make more of these (with specific names and the required parameters) for the individual branches as needed.
		void FillTree_PhotonInformation(G4int input_photontrackid, G4ThreeVector inputvector, G4double input_photonenergy, G4double input_pathlength, G4double input_tracktime);
		void FillTree_BoxInformation(G4double input_x, G4double input_y,G4double input_z);
		void FillTree_ClockTime(string inputclock_time);
		void FillTree_GpsData(G4String inputgps_data);
		void FillTree_MacroPhotonData();


		void CloseTree();

		G4double GetRunID(){return RunID_mng;G4cout << "runID = " << RunID_mng << G4endl; }
		void SetRunID(G4double ID){RunID_mng = ID; }

		G4double GetEventID(){return EventID_mng; G4cout << "eventID = " << EventID_mng << G4endl; }
		void SetEventID(G4double ID){EventID_mng = ID;}

		G4double GetTotalPhotons(){return totalPhotons_mng;}
		void SetTotalPhotons(G4double total){totalPhotons_mng = total;}

		void AddDetectedPhoton(){photonDetectCounter_mng++;}
		void AddAbsorbedPhoton(){photonAbsorbCounter_mng++;}
		void AddEscapedPhoton(){photonEscapeCounter_mng++;}

		void InitialiseCounters(){	photonDetectCounter_mng=0;
									photonAbsorbCounter_mng=0;
									photonEscapeCounter_mng=0;
									overallBounceCounter_mng=0;
									overallEnergyCounter_mng=0;
									overallLengthCounter_mng=0;
									overallLifetimeCounter_mng=0;}

		void AddPhotonBounce(){photonBounceCounter_mng++;}
		void AddBounces(){overallBounceCounter_mng+=photonBounceCounter_mng;}
		
		void InitialiseBounceCounter(){photonBounceCounter_mng=0;}

		void AddEnergyLengthAndTime(G4double energy, G4double length, G4double lifetime){	overallEnergyCounter_mng+=energy;
																							overallLengthCounter_mng+=length;
																							overallLifetimeCounter_mng+=lifetime;}

		void TurnOffBatchInputs(){batchSpecialInputs_mng=false;}
		void TurnOnBatchInputs(){batchSpecialInputs_mng=true;}
		G4bool UsingBatchInputs(){return batchSpecialInputs_mng;}

		void SetBatchGeometry(G4double x,G4double y,G4double z){batchGeomX_mng=x;batchGeomY_mng=y;batchGeomZ_mng=z;}
		G4double GetBatchGeomX(){return batchGeomX_mng;}
		G4double GetBatchGeomY(){return batchGeomY_mng;}
		G4double GetBatchGeomZ(){return batchGeomZ_mng;}


		//Define Get and Set methods for each _mng variable: can be defined inline.

		void PrintToScreen(G4String input) {G4cout << input << G4endl;}



		~WLSRootManager();

};

#endif
