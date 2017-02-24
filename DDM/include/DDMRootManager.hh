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
#include "TGraph.h"
#include "TGraph2D.h"
//#include "TCanvas.h"
#include <sstream>
#include "TFitResult.h"
#include "G4PhysicalConstants.hh"

using namespace std;

class DDMRootManager
{
	private:
		TFile* root_file;
		// Define _mng variables/arrays as needed.
		Double_t TimeStepData_mng[5] = {0};
		Double_t ElectronData_mng[3] = {0};
		Double_t RecoData_mng[4] = {0};
		Int_t EventCounter_mng = -1;
		Int_t ElectronCounter_mng = 0;
		Int_t SecondaryScintCounter_mng = 0;
		Double_t DriftVelocity_mng = 0;
		Double_t ElectricField_mng = 0;
		Double_t Temperature_mng = 0;
		Double_t GasPressure_mng = 0;
		Double_t TankHeight_mng = 0;
		Double_t RecoResults_mng[9] = {0};
		//Double_t InitialMomentum_mng[4] = {0};
		Bool_t IsFirstStep_mng = false;
		Double_t TrueTheta_mng = 0;
		Double_t TruePhi_mng = 0;
		
	public:
		static void CreateRootManager(G4String filename);
		static DDMRootManager* GetRootManager();
		static void DestroyRootManager();

		DDMRootManager(G4String filename);

		void InitialiseTrees();// Add more trees as arguments if desired.
		void InitialiseResultsTree();

		//void FillTree(); Make more of these (with specific names and the required parameters) for the individual branches as needed.
		void FillTree_TimeStepData(G4double input_time, G4double input_x, G4double input_y, G4double input_z, G4double input_energy);
		void FillTree_ElectronData(G4double input_time, G4double input_x, G4double input_y);
		//void FillHist_ElectronGen(G4double input_x, G4double input_y, G4double input_z);
		void FillTree_RecoTrack(G4double input_x, G4double input_y, G4double input_z, G4double input_time);
		//void FillHist_RecoTrack(G4double input_x, G4double input_y, G4double input_z);
		void FillGraph_RecoTrackXY(G4double input_x, G4double input_y);
		void FillGraph_RecoTrackXZ(G4double input_x, G4double input_z);
		void FillGraph_RecoTrackYZ(G4double input_y, G4double input_z);
		void FillGraph_RecoTrack(G4double input_x, G4double input_y, G4double input_z);
		void FillTree_RecoResults(G4double input_tanphi, G4double input_tantheta_xz, G4double input_tantheta_yz);
		void FillGraph_ElectronGen(G4double input_x, G4double input_y, G4double input_z);
	
		void FinaliseEvent();
		void CloseResultsTree();
	
		//void NewBranch();
	
		Double_t CalculateTanThetaFromXZ(Double_t input_tanphi, Double_t input_tanalpha);
		Double_t CalculateTanThetaFromYZ(Double_t input_tanphi, Double_t input_tanbeta);
	
		G4double CalculateSigmaT(G4double input_time, G4double input_mu, G4double input_T);
		G4double CalculateSigmaL(G4double input_time, G4double input_mu, G4double input_T);
	
		G4double CalculateDriftVelocity();

		// define Get and Set methods for each _mng variable: can be defined inline.
		Int_t GetEventCounter() {return EventCounter_mng;}
		
		void SetDriftVelocity(Double_t input) {DriftVelocity_mng = input;}
		Double_t GetDriftVelocity() {return DriftVelocity_mng;}
	
		void SetTankHeight(Double_t input) {TankHeight_mng = input;}
		Double_t GetTankHeight() {return TankHeight_mng;}

		void PrintToScreen(G4String input) {G4cout << input << G4endl;}
	
		// initialise, increment and get counters
		void InitialiseElectronCounter() {ElectronCounter_mng = 0;}
		void IncrementElectronCounter() {ElectronCounter_mng++;}
		Int_t GetElectronCounter() {return ElectronCounter_mng;}
	
		void InitialiseSecondaryScintCounter() {SecondaryScintCounter_mng = 0;}
		void AddToSecondaryScintCounter(Int_t input) {SecondaryScintCounter_mng += input;}
		Int_t GetSecondaryScintCounter_mng() {return SecondaryScintCounter_mng;}
	
		/*void SetInitialMomentum(Double_t input_pabs, Double_t input_px, Double_t input_py, Double_t input_pz)
											{InitialMomentum_mng[0] = input_pabs;
											 InitialMomentum_mng[1] = input_px;
											 InitialMomentum_mng[2] = input_py;
											 InitialMomentum_mng[3] = input_pz;}
		Double_t GetInitialMomentumTanPhi(){return InitialMomentum_mng[2]/InitialMomentum_mng[1];}
		Double_t GetInitialMomentumCosTheta(){return InitialMomentum_mng[3]/InitialMomentum_mng[0];}*/
	
		void SetTrueTheta(Double_t input){TrueTheta_mng = input;}
		void SetTruePhi(Double_t input){TruePhi_mng = input;}
	
		void FlagFirstStep(){IsFirstStep_mng = true;}
		void UnflagFirstStep(){IsFirstStep_mng = false;}
		Bool_t IsFirstStep(){return IsFirstStep_mng;}
	
		void SetElectricField(Double_t input){ElectricField_mng = input;}
		Double_t GetElectricField(){return ElectricField_mng;}
	
		Double_t GetElectronMobility(){return DriftVelocity_mng/ElectricField_mng;}
	
		void SetTemperature(Double_t input){Temperature_mng = input;}
		Double_t GetTemperature(){return Temperature_mng;}
	
		void SetGasPressure(Double_t input){GasPressure_mng = input;}
		Double_t GetGasPressure(){return GasPressure_mng;}

		~DDMRootManager();
};

#endif
