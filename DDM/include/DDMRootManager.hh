#ifndef ROOT_MANAGER_HH
#define ROOT_MANAGER_HH

#include "TFile.h"
#include "TTree.h"
#include "globals.hh"
#include <string>
#include "G4String.hh"
#include "TH1D.h"
#include "TH3I.h"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
//#include "TCanvas.h"
#include <sstream>

using namespace std;

class DDMRootManager
{
	private:
		TFile* root_file;
		//Define _mng variables/arrays as needed.
		Double_t TimeStepData_mng[5] = {0};
		Double_t ElectronData_mng[3] = {0};
		Int_t EventCounter_mng = -1;
		
	public:
		static void CreateRootManager(G4String filename);
		static DDMRootManager* GetRootManager();
		static void DestroyRootManager();

		DDMRootManager(G4String filename);

		void InitialiseTrees();//Add more trees as arguments if desired.

		//void FillTree(); Make more of these (with specific names and the required parameters) for the individual branches as needed.
		void FillTree_TimeStepData(G4double input_time, G4double input_x, G4double input_y, G4double input_z, G4double input_energy);
		void FillTree_ElectronData(G4double input_time, G4double input_x, G4double input_z);
		void FillHist_ElectronGen(G4double input_x, G4double input_y, G4double input_z);
	
		void CloseTrees();
	
		//void NewBranch();

		//Define Get and Set methods for each _mng variable: can be defined inline.
		Int_t GetEventCounter() {return EventCounter_mng;}

		void PrintToScreen(G4String input) {G4cout << input << G4endl;}

		~DDMRootManager();
};

#endif
