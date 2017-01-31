#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* energyLoss_tree;

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

void DDMRootManager::InitialiseTree(G4String treename1)
{
	energyLoss_tree = new TTree("T", treename1.c_str());

	// # BRANCH
	energyLoss_tree -> Branch("energyLoss_data",&TimeStepData_mng, "Time:posx:posy:posz:IonisationEnergy_ev");
	
}

void DDMRootManager::FillTree_TimeStepData(G4int input_time, G4ThreeVector input_vector, G4double input_energy)
{
	TimeStepData_mng[0]=input_time/ms;
	TimeStepData_mng[1]=input_vector.x()/mm;
	TimeStepData_mng[2]=input_vector.y()/mm;
	TimeStepData_mng[3]=input_vector.z()/mm;
	TimeStepData_mng[4]=input_energy/eV;

	energyLoss_tree->Fill(); //Will this fill this branch?

}

void DDMRootManager::CloseTree()
{
	energyLoss_tree->Write();
	
	delete energyLoss_tree;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;

}
