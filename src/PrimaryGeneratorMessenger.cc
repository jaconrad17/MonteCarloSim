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
/// \file optical/OpNovice2/src/PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* Gun)
  : G4UImessenger()
  , fPrimaryAction(Gun)
{
  fGunDir = new G4UIdirectory("/opnovice2/generator/");
  fGunDir->SetGuidance("PrimaryGenerator control");
    
  fSetModeCmd = new G4UIcmdWithAnInteger("/opnovice2/generator/Mode",this);
  fSetModeCmd->SetGuidance("Set the mode of the generator (0 for BEAM or 1 for ROOT");
    
  fSetNEventsCmd = new G4UIcmdWithAnInteger("/opnovice2/generator/Nevents",this);
  fSetNEventsCmd->SetGuidance("Set the numnber of primary events to be generated");
    
  fPolarCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/generator/optPhotonPolar", this);
  fPolarCmd->SetGuidance("Set linear polarization");
  fPolarCmd->SetGuidance("  angle w.r.t. (k,n) plane");
  fPolarCmd->SetParameterName("angle", true);
  fPolarCmd->SetUnitCategory("Angle");
  fPolarCmd->SetDefaultValue(-360.0);
  fPolarCmd->SetDefaultUnit("deg");
  fPolarCmd->AvailableForStates(G4State_Idle);
    
  fSetInputCmd = new G4UIcmdWithAString("/opnovice2/generator/InputFile",this);
  fSetInputCmd->SetGuidance("Set the full name and path of the file with the input ROOT ntuple (in ROOT mode)");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fGunDir;
  delete fSetModeCmd;
  delete fSetNEventsCmd;
  delete fPolarCmd;
  delete fSetInputCmd;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fSetModeCmd)
  {
    fPrimaryAction->SetMode(fSetModeCmd->GetNewIntValue(newValue));
  }
  if(command == fSetNEventsCmd)
  {
    fPrimaryAction->SetNEvents(fSetNEventsCmd->GetNewIntValue(newValue));
  }
  if(command == fSetInputCmd)
  {
    fPrimaryAction->SetUpROOTInput(static_cast<TString>(newValue));
  }
  if(command == fPolarCmd)
  {
    G4double angle = fPolarCmd->GetNewDoubleValue(newValue);
    if(angle == -360.0 * deg)
    {
      fPrimaryAction->SetOptPhotonPolar();
    }
    else
    {
      fPrimaryAction->SetOptPhotonPolar(angle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
