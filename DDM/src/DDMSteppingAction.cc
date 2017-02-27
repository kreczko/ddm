//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DDMSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file DDMSteppingAction.cc
/// \brief Implementation of the DDMSteppingAction class

#include "DDMSteppingAction.hh"
#include "DDMRootManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DDMSteppingAction::DDMSteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DDMSteppingAction::~DDMSteppingAction()
{
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DDMSteppingAction::UserSteppingAction(const G4Step* step)
{
  DDMRootManager* root_manager = DDMRootManager::GetRootManager();
  
  G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     fEventNumber = eventNumber;
     fScintillationCounter = 0;
     fCerenkovCounter = 0;
  }

  G4Track* track = step->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();

  if (ParticleName == "opticalphoton") return;

  const std::vector<const G4Track*>* secondaries =
                                            step->GetSecondaryInCurrentStep();

  if (secondaries->size()>0) {
     for(unsigned int i=0; i<secondaries->size(); ++i) 
     {
        if (secondaries->at(i)->GetParentID()>0)
        {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition())
           {
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Scintillation")
              {
                fScintillationCounter++;
              }
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Cerenkov")
              {
                fCerenkovCounter++;
              }
           }
            else if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4Electron::ElectronDefinition())
            {
              /*if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
               == "Ionisation")*/fIonisationCounter++;
            }
        }
     }
  }
  
  if (track->GetTrackID() == 1)
  {
    // check for first step to log initial momentum
    if (root_manager->IsFirstStep())
    {
      root_manager->UnflagFirstStep();
      G4ThreeVector initialMom = step->GetPreStepPoint()->GetMomentum();
      //root_manager->SetInitialMomentum(initialMom.r(), initialMom.x(), initialMom.y(), initialMom.z());
      root_manager->SetTrueTheta(initialMom.getTheta());
      root_manager->SetTruePhi(initialMom.getPhi());
      //G4cout << "Initial G4ThreeVector theta: " << initialMom.getTheta() << ", phi: " << initialMom.getPhi() << G4endl;
    }
    
    // get time, ionisation energy and position of step
    G4double ionisationEnergy = step->GetTotalEnergyDeposit() - step->GetNonIonizingEnergyDeposit();
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4double stepTime = step->GetPreStepPoint()->GetLocalTime();
    G4cout << "time: " << stepTime << ", pos: " << position << ", E: " << ionisationEnergy << G4endl;
    
    // fill tree
    //root_manager->FillTree_TimeStepData(stepTime,position,ionisationEnergy);
    root_manager->FillTree_TimeStepData(stepTime, position.x(), position.y(), position.z(), ionisationEnergy);
    
    // generate electrons
    G4int numOfElectrons = ionisationEnergy/(15.75962*eV);
    fIonisationCounter += numOfElectrons;
    G4cout << "Electrons manually generated in this step: " << numOfElectrons << G4endl;
    
    // propagate and save electrons
    for (G4int i = 0; i < numOfElectrons; i++)
    {
      root_manager->IncrementElectronCounter();
     
      // *****  electron initial positions  ***********************************************
      
      // randomly generate step fraction
      G4double stepFraction = G4UniformRand();
      
      // set initial x, y and z coords
      G4double initial_x = step->GetPreStepPoint()->GetPosition().x()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().x() - step->GetPreStepPoint()->GetPosition().x());
      G4double initial_y = step->GetPreStepPoint()->GetPosition().y()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().y() - step->GetPreStepPoint()->GetPosition().y());
      G4double initial_z = step->GetPreStepPoint()->GetPosition().z()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().z() - step->GetPreStepPoint()->GetPosition().z());
      G4double initial_time = step->GetPreStepPoint()->GetLocalTime()
                         + stepFraction*(step->GetPostStepPoint()->GetLocalTime() - step->GetPreStepPoint()->GetLocalTime());
      
      //root_manager->FillHist_ElectronGen(initial_x, initial_y, initial_z);
      root_manager->FillGraph_ElectronGen(initial_x, initial_y, initial_z);
      
      // *****  electron drift  ***********************************************************
      
      G4double tankHeight = root_manager->GetTankHeight();
      G4double driftVelocity = root_manager->GetDriftVelocity();
      G4cout << "Drift velocity from manager: " << root_manager->GetDriftVelocity() << G4endl;
      G4double distanceToDrift = tankHeight - initial_z;
      G4double timeToDrift = distanceToDrift/driftVelocity;
      G4double final_time = initial_time + timeToDrift;
      G4double final_x = initial_x;
      G4double final_y = initial_y;
      G4double final_z = tankHeight;
      
      // *****  electron diffusion  ******************************************************
      
      // calculate standard deviation for transerse and longitudinal diffusion
      G4double sigmaT = root_manager->CalculateSigmaT(timeToDrift,
                                                      root_manager->GetElectronMobility(),             
                                                      root_manager->GetTemperature());
      
      G4double sigmaL = root_manager->CalculateSigmaL(timeToDrift, 
                                                      root_manager->GetElectronMobility(),
                                                      root_manager->GetTemperature());
      
      G4cout << "Sigma_T: " << sigmaT << G4endl;
      
      // add Gaussian terms to x and y coords
      final_x += G4RandGauss::shoot(0.0, sigmaT)*mm;
      final_y += G4RandGauss::shoot(0.0, sigmaT)*mm;
      final_z += G4RandGauss::shoot(0.0, sigmaL)*mm;
      
      // ***** scintillation ************************************************************
      
      G4int scintPhotons = root_manager->GetSecondaryScintYield();
      
      G4double scintToCameraDistance = 0.5*m;
      
      for (G4int j=0; j < scintPhotons; j++)
      {
        G4double photonTheta = G4UniformRand()*M_PI;
        G4double photonPhi = G4UniformRand()*2*M_PI;
        
        if (photonTheta < 0.5*M_PI)
        {
          G4double camera_x = final_x + (scintToCameraDistance*tan(photonTheta)*cos(photonPhi));
          G4double camera_y = final_y + (scintToCameraDistance*tan(photonTheta)*sin(photonPhi));
        }
        
        root_manager->FillHist_Camera(final_x, final_y);
      }
      
      // *****  filling data  ************************************************************
      
      // fill data
      root_manager->FillTree_ElectronData(final_time, initial_x, initial_y);  // may need to change time input
      
      // recorded and reconstructed information
      G4double recorded_x = final_x;
      G4double recorded_y = final_y;
      G4double recorded_time = final_time;
      G4double reconstructed_z = tankHeight - (timeToDrift*driftVelocity);
      
      // pixellated view of photons
      root_manager->FillHist_DirectScint(final_x, final_y, scintPhotons);
      
      // fill data to reconstruct track
      //root_manager->FillHist_RecoTrack(recorded_x, recorded_y, reconstructed_z);
      root_manager->FillTree_RecoTrack(recorded_x, recorded_y, reconstructed_z, recorded_time);
      root_manager->FillGraph_RecoTrackXY(recorded_x, recorded_y);
      root_manager->FillGraph_RecoTrack(recorded_x, recorded_y, reconstructed_z);
      root_manager->FillGraph_RecoTrackXZ(recorded_x, reconstructed_z);
      root_manager->FillGraph_RecoTrackYZ(recorded_y, reconstructed_z);
    }
  }
  
  //G4cout << "TrackID: " << track->GetTrackID() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
