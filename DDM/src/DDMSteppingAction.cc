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
      // randomly generate step fraction
      G4double stepFraction = G4UniformRand();
      
      // set initial x, y and z coords
      G4double initial_x = step->GetPreStepPoint()->GetPosition().x()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().x() - step->GetPreStepPoint()->GetPosition().x());
      G4double initial_y = step->GetPreStepPoint()->GetPosition().y()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().y() - step->GetPreStepPoint()->GetPosition().y());
      G4double initial_z = step->GetPreStepPoint()->GetPosition().x()
                         + stepFraction*(step->GetPostStepPoint()->GetPosition().z() - step->GetPreStepPoint()->GetPosition().z());
      G4double initial_time = step->GetPreStepPoint()->GetLocalTime()
                         + stepFraction*(step->GetPostStepPoint()->GetLocalTime() - step->GetPreStepPoint()->GetLocalTime());
      
      root_manager->FillHist_ElectronGen(initial_x, initial_y, initial_z);
      
      // Drift
      G4double tankHeight = root_manager->GetTankHeight();
      G4double driftVelocity = root_manager->GetDriftVelocity();
      G4double distanceToDrift = tankHeight - initial_y;
      G4double timeToDrift = distanceToDrift/driftVelocity;
      G4double final_time = initial_time + timeToDrift;
      G4double final_x = initial_x;
      G4double final_y = tankHeight;
      G4double final_z = initial_z;
      
      // Diffusion
      
      
      // input data
      root_manager->FillTree_ElectronData(final_time, final_x, final_z);
      
      G4double recorded_x = final_x;
      G4double recorded_z = final_z;
      G4double recorded_time = final_time;
      G4double reconstructed_y = tankHeight - (final_time*driftVelocity);

      root_manager->FillHist_RecoTrack(recorded_x, reconstructed_y recorded_z);
      
    }
    
  }
  
  //G4cout << "TrackID: " << track->GetTrackID() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
