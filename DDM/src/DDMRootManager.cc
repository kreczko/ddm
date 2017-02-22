#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* trueTrack_tree;
TTree* electronData_tree;
TTree* recoTrack_tree;
TTree* recoResults_tree;
//TCanvas* c1;
//TH3I* electronGen_hist;
//TH3I* recoTrack_hist;
TGraph* recoTrackXY_graph;
TGraph* recoTrackXZ_graph;
TGraph* recoTrackYZ_graph;
TGraph2D* recoTrack_graph;
TGraph2D* electronGen_graph;

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
	
	// *****  true track information  *************************************
	
	// trueTrackTree
	stringstream trueTrack_treename;
	trueTrack_treename << "trueTrack_" << EventCounter_mng;
	trueTrack_tree = new TTree(trueTrack_treename.str().c_str(), trueTrack_treename.str().c_str());

	trueTrack_tree -> Branch("trueTrack_branch",&TimeStepData_mng, "Time_ns/D:posx_m/D:posy_m/D:posz_m/D:IonisationEnergy_keV/D");
	
	// ***** electron information *****************************************
	
	// electronData_tree
	stringstream electronData_treename;
	electronData_treename << "electronData_" << EventCounter_mng;
	electronData_tree = new TTree(electronData_treename.str().c_str(), electronData_treename.str().c_str());
	
	electronData_tree -> Branch("electronData_branch",&ElectronData_mng, "Time_ns/D:posx_m/D:posz_m/D");
	
	// electronGen_hist
	/*stringstream electronGen_histname;
	electronGen_histname << "electronGen_" << EventCounter_mng;
	electronGen_hist = new TH3I(electronGen_histname.str().c_str(), "Electron generation", 200, -1.0, 1.0, 200, -1.0, 1.0, 200, -1.0, 1.0);
	*/
	// *****  reconstructed track information  *****************************
	
	// recoTrack_tree
	stringstream recoTrack_treename;
	recoTrack_treename << "recoTrack_" << EventCounter_mng;
	recoTrack_tree = new TTree(recoTrack_treename.str().c_str(), recoTrack_treename.str().c_str());
	
	recoTrack_tree -> Branch("recoTrack_branch",&RecoData_mng, "posx_m/D:posy_m/D:posz_m/D:Time_ns/D");
	
	// recoTrack_hist
	/*stringstream recoTrack_histname;
	recoTrack_histname << "recoTrack_H_" << EventCounter_mng;
	recoTrack_hist = new TH3I(recoTrack_histname.str().c_str(), "Reconstructed positions", 200, -1.0, 1.0, 200, -1.0, 1.0, 200, -1.0, 1.0);
	*/
	// recoTrackXY_graph
	recoTrackXY_graph = new TGraph(1);
	stringstream recoTrackXY_graphname;
	recoTrackXY_graphname << "recoTrack_Gxy_" << EventCounter_mng;
	recoTrackXY_graph->SetTitle(recoTrackXY_graphname.str().c_str());
	recoTrackXY_graph->SetName(recoTrackXY_graphname.str().c_str());
	
	// recoTrackXZ_graph
	recoTrackXZ_graph = new TGraph(1);
	stringstream recoTrackXZ_graphname;
	recoTrackXZ_graphname << "recoTrack_Gxz_" << EventCounter_mng;
	recoTrackXZ_graph->SetTitle(recoTrackXZ_graphname.str().c_str());
	recoTrackXZ_graph->SetName(recoTrackXZ_graphname.str().c_str());
	
	// recoTrackYZ_graph
	recoTrackYZ_graph = new TGraph(1);
	stringstream recoTrackYZ_graphname;
	recoTrackYZ_graphname << "recoTrack_Gyz_" << EventCounter_mng;
	recoTrackYZ_graph->SetTitle(recoTrackYZ_graphname.str().c_str());
	recoTrackYZ_graph->SetName(recoTrackYZ_graphname.str().c_str());
	
	// recoTrack_graph (3D graph)
	recoTrack_graph = new TGraph2D(1);
	stringstream recoTrack_graphname;
	recoTrack_graphname << "recoTrack_G_" << EventCounter_mng;
	recoTrack_graph->SetTitle(recoTrack_graphname.str().c_str());
	recoTrack_graph->SetName(recoTrack_graphname.str().c_str());
	
	//electronGen_graph (3D graph)
	electronGen_graph = new TGraph2D(1);
	stringstream electronGen_graphname;
	electronGen_graphname << "electronGen_G_" << EventCounter_mng;
	electronGen_graph->SetTitle(electronGen_graphname.str().c_str());
	electronGen_graph->SetName(electronGen_graphname.str().c_str());
}

void DDMRootManager::InitialiseResultsTree()
{
	recoResults_tree = new TTree("recoResults", "recoResults");
	recoResults_tree -> Branch("recoResults_branch",&RecoResults_mng,
				   "EventNo/D:phi_true/D:theta_true/D:tanphi/D:phi/D:tantheta_xz/D:theta_xz/D:tantheta_yz/D:theta_yz/D");
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

/*void DDMRootManager::FillHist_ElectronGen(G4double input_x, G4double input_y, G4double input_z)
{
	electronGen_hist->Fill(input_x/m, input_y/m, input_z/m);
}*/

void DDMRootManager::FillTree_ElectronData(G4double input_time, G4double input_x, G4double input_y)
{
	ElectronData_mng[0]=input_time/ns;
	ElectronData_mng[1]=input_x/m;
	ElectronData_mng[2]=input_y/m;
	
	electronData_tree->Fill();
}

void DDMRootManager::FillGraph_ElectronGen(G4double input_x, G4double input_y, G4double input_z)
{
	electronGen_graph->Set(ElectronCounter_mng);
	electronGen_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m, input_z/m);
}

/*void DDMRootManager::FillHist_RecoTrack(G4double input_x, G4double input_y, G4double input_z)
{
	recoTrack_hist->Fill(input_x/m, input_y/m, input_z/m);
}*/

void DDMRootManager::FillTree_RecoTrack(G4double input_x, G4double input_y, G4double input_z, G4double input_time)
{
	RecoData_mng[0] = input_x/m;
	RecoData_mng[1] = input_y/m;
	RecoData_mng[2] = input_z/m;
	RecoData_mng[3] = input_time/ns;
	
	recoTrack_tree->Fill();
}

void DDMRootManager::FillGraph_RecoTrackXY(G4double input_x, G4double input_y)
{
	recoTrackXY_graph->Set(ElectronCounter_mng); // Increase size of TGraph as needed
	recoTrackXY_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m); // Set coords of newly-created point
}

void DDMRootManager::FillGraph_RecoTrackYZ(G4double input_y, G4double input_z)
{
	recoTrackYZ_graph->Set(ElectronCounter_mng);
	recoTrackYZ_graph->SetPoint(ElectronCounter_mng - 1, input_y/m, input_z/m);
}

void DDMRootManager::FillGraph_RecoTrackXZ(G4double input_x, G4double input_z)
{
	recoTrackXZ_graph->Set(ElectronCounter_mng);
	recoTrackXZ_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_z/m);
}

void DDMRootManager::FillGraph_RecoTrack(G4double input_x, G4double input_y, G4double input_z)
{
	recoTrack_graph->Set(ElectronCounter_mng);
	recoTrack_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m, input_z/m);
}

void DDMRootManager::FillTree_RecoResults(G4double input_tanphi, G4double input_tantheta_xz, G4double input_tantheta_yz)
{
	RecoResults_mng[0] = EventCounter_mng;
	RecoResults_mng[1] = TruePhi_mng;
	RecoResults_mng[2] = TrueTheta_mng;
	RecoResults_mng[3] = input_tanphi;
	RecoResults_mng[4] = atan(input_tanphi);
	RecoResults_mng[5] = input_tantheta_xz;
	RecoResults_mng[6] = atan(input_tantheta_xz);
	RecoResults_mng[7] = input_tantheta_yz;
	RecoResults_mng[8] = atan(input_tantheta_yz);
	
	recoResults_tree->Fill();
}

G4double CalculateSigmaT(G4double input_time, G4double input_mu, G4double input_T)
{
	G4double transverseD = input_mu*input_T*k_Boltzmann/electron_charge;
	
	G4double sigmaTransverse = sqrt(2*transverseD*input_time);
	return sigmaTransverse;
}

void DDMRootManager::FinaliseEvent()
{
	trueTrack_tree->Write();
	electronData_tree->Write();
	electronGen_graph->Write();
	//electronGen_hist->Write();
	recoTrack_tree->Write();
	//recoTrack_hist->Write();
	
	
	// linear fits
	TFitResultPtr fitXY = recoTrackXY_graph->Fit("pol1", "S");
	TFitResultPtr fitXZ = recoTrackXZ_graph->Fit("pol1", "S");
	TFitResultPtr fitYZ = recoTrackYZ_graph->Fit("pol1", "S");
	
	// fill results tree
	Double_t tanThetaXZ = CalculateTanThetaFromXZ(fitXY->Parameter(1), fitXZ->Parameter(1));
	Double_t tanThetaYZ = CalculateTanThetaFromYZ(fitXY->Parameter(1), fitYZ->Parameter(1));
	
	FillTree_RecoResults(fitXY->Parameter(1), tanThetaXZ, tanThetaYZ);
	
	// label recoTrack_graph axes
	recoTrack_graph->GetXaxis()->SetTitle("x (m)");
	recoTrack_graph->GetYaxis()->SetTitle("y (m)");
	recoTrack_graph->GetZaxis()->SetTitle("z (m)");
	
	// label recoTrackXY_graph axes
	recoTrackXY_graph->GetXaxis()->SetTitle("x (m)");
	recoTrackXY_graph->GetYaxis()->SetTitle("y (m)");
	
	// label recoTrackXZ_graph axes
	recoTrackXZ_graph->GetXaxis()->SetTitle("x (m)");
	recoTrackXZ_graph->GetYaxis()->SetTitle("z (m)");
	
	// label recoTrackYZ_graph axes
	recoTrackYZ_graph->GetXaxis()->SetTitle("y (m)");
	recoTrackYZ_graph->GetYaxis()->SetTitle("z (m)");
	
	recoTrackXY_graph->Write();
	recoTrackXZ_graph->Write();
	recoTrackYZ_graph->Write();
	
	recoTrack_graph->Write();
	
	delete trueTrack_tree;
	delete electronData_tree;
	//delete electronGen_hist;
	delete electronGen_graph;
	delete recoTrack_tree;
	//delete recoTrack_hist;
	delete recoTrackXY_graph;
	delete recoTrackXZ_graph;
	delete recoTrackYZ_graph;
	delete recoTrack_graph;
}

Double_t DDMRootManager::CalculateTanThetaFromXZ(Double_t input_tanphi, Double_t input_tanalpha)
{
	Double_t phi = atan(input_tanphi);
	Double_t tantheta = 1.0/(input_tanalpha * cos(phi));
	return tantheta;
}

Double_t DDMRootManager::CalculateTanThetaFromYZ(Double_t input_tanphi, Double_t input_tanbeta)
{
	Double_t phi = atan(input_tanphi);
	Double_t tantheta = 1.0/(input_tanbeta * sin(phi));
	return tantheta;
}

void DDMRootManager::CloseResultsTree()
{
	recoResults_tree->Write();
	delete recoResults_tree;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;
}
