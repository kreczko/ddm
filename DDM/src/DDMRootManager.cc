#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* trueTrack_tree;
TTree* electronData_tree;
//TCanvas* c1;
TH3I* electronGen_hist;

void DDMRootManager::CreateRootManager(G4String filename)
{
	the_root_manager = new DDMRootManager(filename);
}

DDMRootManager* DDMRootManager::GetRootManager()
{
	return the_root_manager;
}

void DDMRootManager::DestroyRootManager()
{
	delete the_root_manager;
}

DDMRootManager::DDMRootManager(G4String filename)
{
	root_file = new TFile(filename.c_str(),"RECREATE");
}

void DDMRootManager::InitialiseTrees()
{
	EventCounter_mng++;
	
	//c1 = new TCanvas("c1");
	
	// True track information:
	stringstream trueTrack_treename;
	trueTrack_treename << "trueTrack_" << EventCounter_mng;
	trueTrack_tree = new TTree(trueTrack_treename.str().c_str(), trueTrack_treename.str().c_str());

	trueTrack_tree -> Branch("trueTrack_branch",&TimeStepData_mng, "Time_ns/D:posx_m/D:posy_m/D:posz_m/D:IonisationEnergy_keV/D");
	
	// Electron information:
	stringstream electronData_treename;
	electronData_treename << "electronData_" << EventCounter_mng;
	electronData_tree = new TTree(electronData_treename.str().c_str(), electronData_treename.str().c_str());
	
	electronData_tree -> Branch("electronData_branch",&ElectronData_mng, "Time_ns/D:posx_m/D:posz_m/D");
	electronGen_hist = new TH3I("electronGen_hist", "Electron generation", 200, -1.0, 1.0, 200, -1.0, 1.0, 200, -1.0, 1.0);
}

void DDMRootManager::FillTree_TimeStepData(G4double input_time, G4double input_x, G4double input_y, G4double input_z, G4double input_energy)
{
	TimeStepData_mng[0]=input_time/ns;
	TimeStepData_mng[1]=input_x/m;
	TimeStepData_mng[2]=input_y/m;
	TimeStepData_mng[3]=input_z/m;
	TimeStepData_mng[4]=input_energy/keV;

	trueTrack_tree->Fill();
}

void DDMRootManager::FillHist_ElectronGen(G4double input_x, G4double input_y, G4double input_z)
{
	electronGen_hist->Fill(input_x/m, input_y/m, input_z/m);
}

void DDMRootManager::FillTree_ElectronData(G4double input_time, G4double input_x, G4double input_z)
{
	ElectronData_mng[0]=input_time/ns;
	ElectronData_mng[1]=input_x/m;
	ElectronData_mng[2]=input_z/m;
	
	electronData_tree->Fill();
}

void DDMRootManager::CloseTrees()
{
	trueTrack_tree->Write();
	electronData_tree->Write();
	electronGen_hist->Write();
	
	
	delete trueTrack_tree;
	delete electronData_tree;
	delete electronGen_hist;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;
}
