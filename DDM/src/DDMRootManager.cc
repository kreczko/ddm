#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* photon_tree;
TTree* meta_tree;
TTree* macro_photon_tree;
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

void DDMRootManager::InitialiseTree(G4String treename1, G4String treename2, G4String treename3)
{
	photon_tree = new TTree("T", treename1.c_str());
	meta_tree = new TTree ("T2", treename2.c_str());
	macro_photon_tree = new TTree ("T3", treename3.c_str());

	// # BRANCH
	photon_tree -> Branch("photon_data",&PhotonInformation_mng, "EventID/D:TrackID:posx:posy:posz:bounces/D:E_ev:pathlength_mm:tracktime_ns");
	meta_tree -> Branch ("box_data", &BoxInformation_mng,"RunID/D:x_dim:y_dim:z_dim");//dimensions of box
	meta_tree -> Branch ("clock_time", &ClockTime_mng,"clock_time");
	meta_tree -> Branch ("gps_data", &GpsData_mng,"gps_data");
	macro_photon_tree -> Branch ("macro_photon_data", &MacroPhotonData_mng,"RunID/D:no_detected/D:no_absorbed_in_cube/D:no_escaped/D:total_num/D:av_no_bounces:av_E_ev:av_pathlength_mm:av_tracktime_ns");

}

void DDMRootManager::FillTree_PhotonInformation(G4int input_photontrackid, G4ThreeVector inputvector, G4double input_photonenergy, G4double input_pathlength, G4double input_tracktime)

{
	PhotonInformation_mng[0]=EventID_mng;
	PhotonInformation_mng[1]=input_photontrackid;
	PhotonInformation_mng[2]=inputvector.x()/mm;
	PhotonInformation_mng[3]=inputvector.y()/mm;
	PhotonInformation_mng[4]=inputvector.z()/mm;
	PhotonInformation_mng[5]=photonBounceCounter_mng;
	PhotonInformation_mng[6]=input_photonenergy/eV;
	PhotonInformation_mng[7]=input_pathlength/mm;
	PhotonInformation_mng[8]=input_tracktime/ns;


	photon_tree->Fill(); //Will this fill this branch?



}

void DDMRootManager::FillTree_BoxInformation(G4double input_x, G4double input_y,G4double input_z)

{

	BoxInformation_mng[0]=RunID_mng;
	BoxInformation_mng[1]=input_x/mm;
	BoxInformation_mng[2]=input_y/mm;
	BoxInformation_mng[3]=input_z/mm;
	


	meta_tree-> Fill();
}

void DDMRootManager::FillTree_ClockTime(string inputclock_time)
{
	ClockTime_mng=inputclock_time;

	meta_tree-> Fill();
}		


void DDMRootManager::FillTree_GpsData(G4String inputgps_data)
{
	GpsData_mng=inputgps_data;

	meta_tree-> Fill();
}

void DDMRootManager::FillTree_MacroPhotonData()

{
	MacroPhotonData_mng[0]=RunID_mng;
	MacroPhotonData_mng[1]=photonDetectCounter_mng;
	MacroPhotonData_mng[2]=photonAbsorbCounter_mng;
	MacroPhotonData_mng[3]=photonEscapeCounter_mng;
	MacroPhotonData_mng[4]=totalPhotons_mng;
	MacroPhotonData_mng[5]=(G4double)overallBounceCounter_mng/(G4double)photonDetectCounter_mng;
	MacroPhotonData_mng[6]=(overallEnergyCounter_mng/(G4double)photonDetectCounter_mng)/eV;
	MacroPhotonData_mng[7]=(overallLengthCounter_mng/(G4double)photonDetectCounter_mng)/mm;
	MacroPhotonData_mng[8]=(overallLifetimeCounter_mng/(G4double)photonDetectCounter_mng)/ns;

	macro_photon_tree->Fill();
}




void DDMRootManager::CloseTree()
{
	photon_tree->Write();
	meta_tree->Write();
	macro_photon_tree -> Write();
	delete photon_tree;
	delete meta_tree;
	delete macro_photon_tree;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;

}
