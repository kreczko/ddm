#include "DDMRootManager.hh"

DDMRootManager* the_root_manager;
TTree* trueTrack_tree;
TTree* electronData_tree;
TTree* recoTrack_tree;
TTree* recoResults_tree;
TTree* recoResultsCamera_tree;
//TCanvas* c1;
//TH3I* electronGen_hist;
//TH3I* recoTrack_hist;
TGraph* recoTrackXY_graph;
TGraph* recoTrackXZ_graph;
TGraph* recoTrackYZ_graph;
TGraph2D* recoTrack_graph;
TGraph2D* electronGen_graph;
TH2I* directScint_hist;
TH2I* camera_hist;
TLine* cameraFitLine;
TH2I* cameraXZ_hist;
TH2I* cameraYZ_hist;

TGraph* fitCamera_graph;

//TH1D* cameraProjectionX_hist;
//TH1D* cameraProjectionY_hist;

void DDMRootManager::CreateRootManager()
{
	the_root_manager = new DDMRootManager();
}

DDMRootManager* DDMRootManager::GetRootManager()
{
	return the_root_manager;
}

void DDMRootManager::DestroyRootManager()
{
	delete the_root_manager;
}

DDMRootManager::DDMRootManager()
{
	/*if (filename == "[dynamic]")
	{
		stringstream filename_stream;
		filename_stream << "ddm_N" << EventCounter_mng << "_p" << GasPressure_mng << "_E" << ElectricField_mng << ".root";
		filename = filename_stream.str();
	}
	root_file = new TFile(filename.c_str(),"RECREATE");*/
}

void DDMRootManager::CreateOutputFile(G4String filename)
{
	if (filename == "[dynamic]")
	{
		stringstream filename_stream;
		filename_stream << "ddm_p" << GasPressure_mng/atmosphere << ".root";
		filename = filename_stream.str();
	}
	root_file = new TFile(filename.c_str(),"RECREATE");
}

void DDMRootManager::ReadParameterFile(G4String filename)
{
	ifstream parameterFile;
	parameterFile.open(filename);
	
	G4String word;
	Int_t paramCount = 0;
	Int_t expectedParamCount = 2;
	
	while (parameterFile >> word)
	{
		if (word == "Pressure:") 
		{
			paramCount++;
			parameterFile >> word;
			SetGasPressure(atof(word)*atmosphere);
			
			G4cout << "Input pressure: " << word << " atm" << G4endl;
		}
		
		if (word == "TimeRes:") 
		{
			paramCount++;
			parameterFile >> word;
			TimeResolution_mng = atof(word)*ns;
			
			G4cout << "Input time resolution: " << word << " ns" << G4endl;
		}
	}
	
	if (paramCount != expectedParamCount)
	{
		G4cout << "WARNING: Parameter file is missing " << expectedParamCount - paramCount << " values!!!" << G4endl;
	}
}

void DDMRootManager::InitialiseTrees()
{
	EventCounter_mng++;
	
	//c1 = new TCanvas("c1");
	
	// ***********************************  true track information  *************************************
	
	// trueTrackTree
	stringstream trueTrack_treename;
	trueTrack_treename << "trueTrack_" << EventCounter_mng;
	trueTrack_tree = new TTree(trueTrack_treename.str().c_str(), trueTrack_treename.str().c_str());

	trueTrack_tree -> Branch("trueTrack_branch",&TimeStepData_mng, "Time_ns/D:posx_m/D:posy_m/D:posz_m/D:IonisationEnergy_keV/D");
	
	// *********************************** electron information *****************************************
	
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
	
	// ***********************  reconstructed track information (perfect detection)  ********************
	
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
	
	// *********************************** scintillation information ************************************
	
	// directScint_hist
	stringstream directScint_histname;
	directScint_histname << "directScint_" << EventCounter_mng;
	directScint_hist = new TH2I(directScint_histname.str().c_str(), "Scintillation photons", 200, -1.0, 1.0, 200, -1.0, 1.0);
	
	// camera_hist
	stringstream camera_histname;
	camera_histname << "camera_" << EventCounter_mng;
	camera_hist = new TH2I(camera_histname.str().c_str(), "Camera image",
			       CameraResolution_mng, -SensorEffectiveX_mng/m, SensorEffectiveX_mng/m,
			       CameraResolution_mng, -SensorEffectiveY_mng/m, SensorEffectiveY_mng/m);
	//camera_hist = new TH2I(camera_histname.str().c_str(), "Camera image", 1000, -1.0, 1.0, 1000, -1.0, 1.0);
	
	// calcualate z range for camera projection histograms
	Int_t SnapshotNumber = (2.0*TankHeight_mng/(TimeResolution_mng*DriftVelocity_mng));
	G4cout << "Calculated snapshot number: " << SnapshotNumber << G4endl;
	Double_t CameraRange_z = SnapshotNumber * TimeResolution_mng * DriftVelocity_mng;
	
	// cameraXZ_hist
	stringstream cameraXZ_histname;
	cameraXZ_histname << "cameraXZ_" << EventCounter_mng;
	cameraXZ_hist = new TH2I(cameraXZ_histname.str().c_str(), "Camera image XZ",
			       CameraResolution_mng, -SensorEffectiveX_mng/m, SensorEffectiveX_mng/m,
			       SnapshotNumber, -0.5*CameraRange_z/m, 0.5*CameraRange_z/m);
	
	// cameraYZ_hist
	stringstream cameraYZ_histname;
	cameraYZ_histname << "cameraYZ_" << EventCounter_mng;
	cameraYZ_hist = new TH2I(cameraYZ_histname.str().c_str(), "Camera image YZ",
			       CameraResolution_mng, -SensorEffectiveY_mng/m, SensorEffectiveY_mng/m,
			       SnapshotNumber, -0.5*CameraRange_z/m, 0.5*CameraRange_z/m);
	
	/*
	// cameraProjectionX_hist
	cameraProjectionX_hist = new TH1D("ProjectionXName", "ProjectionX",
					  CameraResolution_mng, -SensorEffectiveX_mng/m, SensorEffectiveX_mng/m);
	
	// cameraProjectionY_hist
	cameraProjectionY_hist = new TH1D("ProjectionYName", "ProjectionY",
					  CameraResolution_mng, -SensorEffectiveY_mng/m, SensorEffectiveY_mng/m);
	*/
}

void DDMRootManager::InitialiseResultsTree()
{
	// old version of results tree
	recoResults_tree = new TTree("recoResults", "recoResults");
	//recoResults_tree -> Branch("recoResults_branch",&RecoResults_mng,
				   //"EventNo/D:phi_true/D:theta_true/D:tanphi/D:phi/D:tantheta_xz/D:theta_xz/D:tantheta_yz/D:theta_yz/D");
	
	recoResults_tree -> Branch("recoResults_branch",&RecoResults_mng,
				   "EventNo/D:phi_true/D:theta_true/D:tanphi/D:phi/D:tantheta_xz/D:theta_xz/D:tantheta_yz/D:theta_yz/D:tanphi_scint/D:phi_scint/D");

	// camera results tree
	recoResultsCamera_tree = new TTree("recoResultsCamera", "recoResultsCamera");
	
	recoResultsCamera_tree -> Branch("recoResultsCamera_branch",&RecoResultsCamera_mng,
				   "EventNo/D:phi_true/D:theta_true/D:camera_phi/D:camera_theta_xz/D:camera_theta_yz/D");
}

void DDMRootManager::FillTree_TimeStepData(Double_t input_time, Double_t input_x, Double_t input_y, Double_t input_z, Double_t input_energy)
{
	TimeStepData_mng[0]=input_time/ns;
	TimeStepData_mng[1]=input_x/m;
	TimeStepData_mng[2]=input_y/m;
	TimeStepData_mng[3]=input_z/m;
	TimeStepData_mng[4]=input_energy/keV;

	trueTrack_tree->Fill();
}

/*void DDMRootManager::FillHist_ElectronGen(Double_t input_x, Double_t input_y, Double_t input_z)
{
	electronGen_hist->Fill(input_x/m, input_y/m, input_z/m);
}*/

void DDMRootManager::FillTree_ElectronData(Double_t input_time, Double_t input_x, Double_t input_y)
{
	ElectronData_mng[0]=input_time/ns;
	ElectronData_mng[1]=input_x/m;
	ElectronData_mng[2]=input_y/m;
	
	electronData_tree->Fill();
}

void DDMRootManager::FillGraph_ElectronGen(Double_t input_x, Double_t input_y, Double_t input_z)
{
	electronGen_graph->Set(ElectronCounter_mng);
	electronGen_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m, input_z/m);
}

/*void DDMRootManager::FillHist_RecoTrack(Double_t input_x, Double_t input_y, Double_t input_z)
{
	recoTrack_hist->Fill(input_x/m, input_y/m, input_z/m);
}*/

void DDMRootManager::FillTree_RecoTrack(Double_t input_x, Double_t input_y, Double_t input_z, Double_t input_time)
{
	RecoData_mng[0] = input_x/m;
	RecoData_mng[1] = input_y/m;
	RecoData_mng[2] = input_z/m;
	RecoData_mng[3] = input_time/ns;
	
	recoTrack_tree->Fill();
}

void DDMRootManager::FillGraph_RecoTrackXY(Double_t input_x, Double_t input_y)
{
	recoTrackXY_graph->Set(ElectronCounter_mng); // Increase size of TGraph as needed
	recoTrackXY_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m); // Set coords of newly-created point
}

void DDMRootManager::FillGraph_RecoTrackYZ(Double_t input_y, Double_t input_z)
{
	recoTrackYZ_graph->Set(ElectronCounter_mng);
	recoTrackYZ_graph->SetPoint(ElectronCounter_mng - 1, input_y/m, input_z/m);
}

void DDMRootManager::FillGraph_RecoTrackXZ(Double_t input_x, Double_t input_z)
{
	recoTrackXZ_graph->Set(ElectronCounter_mng);
	recoTrackXZ_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_z/m);
}

void DDMRootManager::FillGraph_RecoTrack(Double_t input_x, Double_t input_y, Double_t input_z)
{
	recoTrack_graph->Set(ElectronCounter_mng);
	recoTrack_graph->SetPoint(ElectronCounter_mng - 1, input_x/m, input_y/m, input_z/m);
}

void DDMRootManager::FillTree_RecoResults(Double_t input_tanphi, Double_t input_tantheta_xz, Double_t input_tantheta_yz, Double_t input_tanphi_scint)
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
	RecoResults_mng[9] = input_tanphi_scint;
	RecoResults_mng[10] = atan(input_tanphi_scint);
	
	recoResults_tree->Fill();
}

void DDMRootManager::FillTree_RecoResultsCamera(Double_t input_tanphi, Double_t input_tantheta_xz, Double_t input_tantheta_yz)
{
	RecoResultsCamera_mng[0] = EventCounter_mng;
	RecoResultsCamera_mng[1] = TruePhi_mng;
	RecoResultsCamera_mng[2] = TrueTheta_mng;
	RecoResultsCamera_mng[3] = atan(input_tanphi);
	RecoResultsCamera_mng[4] = atan(input_tantheta_xz);
	RecoResultsCamera_mng[5] = atan(input_tantheta_yz);
	
	recoResultsCamera_tree->Fill();
}

void DDMRootManager::FillHist_DirectScint(Double_t input_x, Double_t input_y)
{
	//directScint_hist->Fill(input_x/m, input_y/m, input_photonNum);
	directScint_hist->Fill(input_x/m, input_y/m);
}

void DDMRootManager::FillHist_Camera(Double_t input_x, Double_t input_y)
{
	camera_hist->Fill(input_x/m, input_y/m);
}

void DDMRootManager::FillHist_CameraXZ(Double_t input_x, Double_t input_z)
{
	cameraXZ_hist->Fill(input_x/m, input_z/m);
}

void DDMRootManager::FillHist_CameraYZ(Double_t input_y, Double_t input_z)
{
	cameraYZ_hist->Fill(input_y/m, input_z/m);
}

Double_t DDMRootManager::FitCameraHist(TH2I* input_hist)
{
	
	//TGraph* fitCamera_graph = new TGraph(1);
	fitCamera_graph = new TGraph(1);
	
	fitCamera_graph->Set(SecondaryScintCounter_mng);
	
	//G4int photonCut = camera_hist->GetMaximum()/4;
	
	G4int point = 0;
	
	G4cout << "Camera histogram size: " << input_hist->GetNbinsX() << "x" << input_hist->GetNbinsY() << G4endl;
	
	input_hist->ClearUnderflowAndOverflow();
	
	// loop through all bins
	for(G4int binx = 1; binx <= input_hist->GetNbinsX(); binx++)
	{
		for(G4int biny = 1; biny <= input_hist->GetNbinsY(); biny++)
		{
			// find bin centres
			Double_t binCentreX = input_hist->GetXaxis()->GetBinCenter(binx);
			Double_t binCentreY = input_hist->GetYaxis()->GetBinCenter(biny);
			
			// check number of entries in bin is more than the cut
			//if (camera_hist->GetBinContent(binx, biny) > photonCut)
			//{
				// plot a point at bin centre for each photon in bin
				for(Int_t photonPerBin = 0; photonPerBin < input_hist->GetBinContent(binx, biny); photonPerBin++)
				{	
					fitCamera_graph->SetPoint(point, binCentreX, binCentreY);
					point++;
				}
			//}
		}
	}

	// linear fit of graph
	TFitResultPtr cameraFit = fitCamera_graph->Fit("pol1", "S");
	
	// put fit gradient into ROOT manager variable
	CameraTanPhi_mng = cameraFit->Parameter(1);

	// calculate start point for fit line
	Double_t start_x = -1.0;
	Double_t start_y = (start_x*cameraFit->Parameter(1)) + cameraFit->Parameter(0);
	G4cout << "Starting point of fit line: " << start_x << ", " << start_y << G4endl;
	
	// calculate end point for fit line
	Double_t end_x = 1.0;
	Double_t end_y = (end_x*cameraFit->Parameter(1)) + cameraFit->Parameter(0);
	G4cout << "Ending point of fit line: " << end_x << ", " << end_y << G4endl;
	
	// overlay fit line onto camera histogram
	cameraFitLine = new TLine(start_x, start_y, end_x, end_y);
	cameraFitLine->SetLineColor(kRed);
	input_hist->GetListOfFunctions()->Add(cameraFitLine);
	PrintToScreen("Fit line drawn onto camera image.");
	
	// save gradient of fit in a variable
	Double_t cameraGradient = cameraFit->Parameter(1);

	// write sacrifical graph to file
	fitCamera_graph->Write();
	
	delete fitCamera_graph;
	
	//G4cout << "Fit of camera image complete." << G4endl;
	
	// return gradient
	return cameraGradient;
}

G4double DDMRootManager::CalculateDriftVelocity()
{
	G4cout << "Calculating drift velocity..." << G4endl;
  	G4double reducedField = ElectricField_mng/GasPressure_mng;
  	G4double reducedField_VperCmTorr = reducedField/(volt*760/(cm*atmosphere));
  	DriftVelocity_mng = pow(reducedField_VperCmTorr, 0.85)*3.0e5*cm/s;
  	G4cout << "Drift velocity = " << DriftVelocity_mng/(cm/s) << G4endl;
	
  	if ((reducedField_VperCmTorr < 2) || (reducedField_VperCmTorr > 1000))
  	{
    		G4cout << "WARNING: reduced field of " << reducedField_VperCmTorr << " V/(cm Torr) is outside reliable linear range. Calculated drift velocity may be inaccurate." << G4endl;
  	}
	
	return DriftVelocity_mng;
}

G4double DDMRootManager::CalculateSecondaryScintYield(Double_t input_avalancheField)
{
	// Y/p (photons electron^-1 cm ^-1 bar^-1) = 81 E/p (kV cm^-1 bar^-1) - 47
	
	G4double gasPressureBar = GasPressure_mng/bar;
	G4double avalancheFieldPerCm = input_avalancheField/(kilovolt/cm);
	
	SecondaryScintYield_mng = gasPressureBar * (81.0 * (avalancheFieldPerCm/gasPressureBar) - 47.0);
	
	if (avalancheFieldPerCm/gasPressureBar < 0.7 || avalancheFieldPerCm/gasPressureBar > 3.0)
	{
		G4cout << "WARNING: reduced field of " << avalancheFieldPerCm/gasPressureBar << " (kV/cm) out of linear range. Scintillation yield NOT reliable." << G4endl;
	}
	
	return SecondaryScintYield_mng;
}

G4double DDMRootManager::CalculateSigmaT(Double_t input_time, Double_t input_mu, Double_t input_T)
{
	G4double transverseD = input_mu*input_T*k_Boltzmann/(-electron_charge);
	
	G4double sigmaTransverse = sqrt(2*transverseD*input_time);
	return sigmaTransverse;
}

G4double DDMRootManager::CalculateSigmaL(Double_t input_time, Double_t input_mu, Double_t input_T)
{
	G4double longitudinalD = input_mu*input_T*k_Boltzmann/(-electron_charge);
	
	G4double sigmaLongitudinal = sqrt(2*longitudinalD*input_time);
	return sigmaLongitudinal;
}

void DDMRootManager::FinaliseEvent()
{
	if (IsStreamliningOff())
	{
		trueTrack_tree->Write();
		electronData_tree->Write();
		electronGen_graph->Write();
		//electronGen_hist->Write();
		recoTrack_tree->Write();
		//recoTrack_hist->Write();
		directScint_hist->Write();
	}
	
	// TLine test
	/*TLine* testLine = new TLine(-0.5, 0.5, 0.5, -0.5);
	testLine->SetLineColor(kRed);
	camera_hist->GetListOfFunctions()->Add(testLine);*/	
	
	if (IsStreamliningOff())
	{
		// linear fits of each track projection
		TFitResultPtr fitXY = recoTrackXY_graph->Fit("pol1", "S");
		TFitResultPtr fitXZ = recoTrackXZ_graph->Fit("pol1", "S");
		TFitResultPtr fitYZ = recoTrackYZ_graph->Fit("pol1", "S");
	
		// calculate tan(theta) from the fits of XZ and YZ projections
		Double_t tanThetaXZ = CalculateTanThetaFromXZ(fitXY->Parameter(1), fitXZ->Parameter(1));
		Double_t tanThetaYZ = CalculateTanThetaFromYZ(fitXY->Parameter(1), fitYZ->Parameter(1));
	
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
		
		FillTree_RecoResults(fitXY->Parameter(1), tanThetaXZ, tanThetaYZ, CameraTanPhi_mng);
	}
	
	// calculate phi from linear fit of camera histogram
	Double_t cameraTanPhi =  FitCameraHist(camera_hist);
	
	// calculate alpha and beta angles from the fits of XZ and YZ camera projections
	Double_t cameraTanAlpha = FitCameraHist(cameraXZ_hist);
	Double_t cameraTanBeta = FitCameraHist(cameraYZ_hist);
	
	// calculate tan(theta) from camera XZ and YZ projections
	Double_t cameraTanTheta_xz = CalculateTanThetaFromXZ(cameraTanPhi, cameraTanAlpha);
	Double_t cameraTanTheta_yz = CalculateTanThetaFromYZ(cameraTanPhi, cameraTanBeta);
	
	//fitCamera_graph->Write();
	
	camera_hist->Write();
	cameraXZ_hist->Write();
	cameraYZ_hist->Write();
	
	// print skewness along x and skewness along y of camera image
	G4cout << "skewness x = " << camera_hist->GetSkewness(1) << G4endl;
	G4cout << "skewness y = " << camera_hist->GetSkewness(2) << G4endl;
	
	Double_t deviation = CalculateVectorAngle(cameraTanPhi, cameraTanTheta_xz, camera_hist->GetSkewness(1));
	G4cout << "Directional deviation: " << deviation << G4endl;
	
	// fill camera results tree
	FillTree_RecoResultsCamera(cameraTanPhi, cameraTanTheta_xz, cameraTanTheta_yz);
	
	/*
	cameraProjectionX_hist = camera_hist->ProjectionX("Camera_projection_X", 1, CameraResolution_mng);
	cameraProjectionY_hist = camera_hist->ProjectionY("Camera_projection_Y", 1, CameraResolution_mng);
	
	cameraProjectionX_hist->Write();
	cameraProjectionY_hist->Write();
	*/
	
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
	delete directScint_hist;
	delete camera_hist;
	delete cameraXZ_hist;
	delete cameraYZ_hist;
	
	//delete fitCamera_graph;
	
	//delete cameraProjectionX_hist;
	//delete cameraProjectionY_hist;
	
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

Double_t DDMRootManager::CalculateVectorAngle(Double_t input_tanphi, Double_t input_tantheta, Double_t input_skewnessX)
{
	G4ThreeVector* trueDirection = new G4ThreeVector(1.0, 0.0, 0.0);
	trueDirection->setPhi(TruePhi_mng);
	trueDirection->setTheta(TrueTheta_mng);
	
	Double_t theta = atan(input_tantheta);
	
	Double_t phi = atan(input_tanphi);
	if (input_skewnessX > 0) {phi += M_PI*rad;} // apply head-tailing, assuming more electrons at head
	
	G4ThreeVector* recoDirection = new G4ThreeVector(1.0, 0.0, 0.0);
	recoDirection->setPhi(phi);
	recoDirection->setTheta(theta);
	
	Double_t deviation = recoDirection->angle(*trueDirection);
	
	return deviation;
}

void DDMRootManager::CloseResultsTree()
{
	if (IsStreamliningOff()) {recoResults_tree->Write();}
	recoResultsCamera_tree->Write();
	delete recoResults_tree;
	delete recoResultsCamera_tree;
}

DDMRootManager::~DDMRootManager()
{
	root_file->Write();
	delete root_file;
}
