#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* trueTrack_tree;
TTree* electronData_tree;
TTree* recoTrack_tree;
//TCanvas* c1;
TH3I* electronGen_hist;
TH3I* recoTrack_hist;

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
	
	// true track information
	stringstream trueTrack_treename;
	trueTrack_treename << "trueTrack_" << EventCounter_mng;
	trueTrack_tree = new TTree(trueTrack_treename.str().c_str(), trueTrack_treename.str().c_str());

	trueTrack_tree -> Branch("trueTrack_branch",&TimeStepData_mng, "Time_ns/D:posx_m/D:posy_m/D:posz_m/D:IonisationEnergy_keV/D");
	
	// electron information
	stringstream electronData_treename;
	electronData_treename << "electronData_" << EventCounter_mng;
	electronData_tree = new TTree(electronData_treename.str().c_str(), electronData_treename.str().c_str());
	
	electronData_tree -> Branch("electronData_branch",&ElectronData_mng, "Time_ns/D:posx_m/D:posz_m/D");
	
	stringstream electronGen_histname;
	electronGen_histname << "electronGen_" << EventCounter_mng;
	electronGen_hist = new TH3I(electronGen_histname.str().c_str(), "Electron generation", 200, -1.0, 1.0, 200, -1.0, 1.0, 200, -1.0, 1.0);
	
	// reconstructed track information
	stringstream recoTrack_treename;
	recoTrack_treename << "recoTrack_" << EventCounter_mng;
	recoTrack_tree = new TTree(recoTrack_treename.str().c_str(), recoTrack_treename.str().c_str());
	
	recoTrack_tree -> Branch("recoTrack_branch",&RecoData_mng, "posx_m/D:posy_m/D:posz_m/D");
	
	stringstream recoTrack_histname;
	recoTrack_histname << "recoTrack_H_" << EventCounter_mng;
	recoTrack_hist = new TH3I(recoTrack_histname.str().c_str(), "Reconstructed positions", 200, -1.0, 1.0, 200, -1.0, 1.0, 200, -1.0, 1.0);
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

void DDMRootManager::FillHist_RecoTrack(G4double input_x, G4double input_y, G4double input_z)
{
	recoTrack_hist->Fill(input_x/m, input_y/m, input_z/m);
}

void DDMRootManager::FillTree_RecoTrack(G4double input_x, G4double input_y, G4double input_z)
{
	RecoData_mng[0] = input_x/m;
	RecoData_mng[1] = input_y/m;
	RecoData_mng[2] = input_z/m;
	
	recoTrack_tree->Fill();
}

void DDMRootManager::CloseTrees()
{
	trueTrack_tree->Write();
	electronData_tree->Write();
	electronGen_hist->Write();
	recoTrack_tree->Write();
	recoTrack_hist->Write();
	
	delete trueTrack_tree;
	delete electronData_tree;
	delete electronGen_hist;
	delete recoTrack_tree;
	delete recoTrack_hist;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;
}
