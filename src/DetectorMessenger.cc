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
/// \file optical/OpNovice2/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4OpticalSurface.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include <sstream>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Detect)
  : G4UImessenger()
  , fDetector(Detect)
{
    fOpticalDir = new G4UIdirectory("/opnovice2/detector/");
    fOpticalDir->SetGuidance("Parameters for optical simulation.");

    fNPSAngleCmd = new G4UIcmdWithADoubleUnit("/opnovice2/detector/setEarmAngle",this);
    fNPSAngleCmd->SetGuidance("Set Earm angle (value and unit");
    fNPSAngleCmd->SetUnitCategory("Angle");
    
    fNPSDistanceCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setEarmDistance",this);
    fNPSDistanceCmd->SetGuidance("Set Earm distance (value and unit)");
    fNPSDistanceCmd->SetUnitCategory("Length");
    
    fNPSShieldThicknessCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setEarmShieldThickness",this);
    fNPSShieldThicknessCmd->SetGuidance("Set Earm lead shield thickness (value and unit)");
    fNPSShieldThicknessCmd->SetUnitCategory("Length");
    
    fHCALAngleCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setHarmAngle",this);
    fHCALAngleCmd->SetGuidance("Set Harm angle (value and unit)");
    fHCALAngleCmd->SetUnitCategory("Angle");
    
    fHCALDistanceCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setHarmDistance",this);
    fHCALDistanceCmd->SetGuidance("Set Harm distance (value and unit");
    fHCALDistanceCmd->SetUnitCategory("Length");
    
    fHCALShieldThicknessCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/SetHarmShieldThickness",this);
    fHCALShieldThicknessCmd->SetGuidance("Set Harm lead shield thickness (value and unit)");
    fHCALShieldThicknessCmd->SetUnitCategory("Length");
    
    fWindowThicknessCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setWindowThickness",this);
    fWindowThicknessCmd->SetGuidance("Set scattering chamber window thickness (value and unit)");
    fWindowThicknessCmd->SetUnitCategory("Length");
    
    fTargetLengthCmd = new G4UIcmdWithADoubleAndUnit("/opnovice2/detector/setTargetLength",this);
    fTargetLengthCmd->SetGuidance("Set LH2 target length (value and unit)");
    fTargetLengthCmd->SetUnitCategory("Length");
    
    fBeamlineCmd = new G4UIcmdWithAnInteger("/opnovice2/detector/setBeamline",this);
    fBeamlineCmd->SetGuidance("Set an integer flag for whether to include scattering chamber and beamline (0 or 1)");
    
    fUpdateCmd = new G4UIcommand("/opnovice2/detector/update",this);
    fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
    fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.")

    fWorldMatPropVectorCmd =
    new G4UIcmdWithAString("/opnovice2/detector/worldProperty", this);
    fWorldMatPropVectorCmd->SetGuidance("Set material property vector ");
    fWorldMatPropVectorCmd->SetGuidance("for the world.");
    fWorldMatPropVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fWorldMatPropVectorCmd->SetToBeBroadcasted(false);

    fWorldMatPropConstCmd =
    new G4UIcmdWithAString("/opnovice2/detector/worldConstProperty", this);
    fWorldMatPropConstCmd->SetGuidance("Set material constant property");
    fWorldMatPropConstCmd->SetGuidance(" for the world.");
    fWorldMatPropConstCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fWorldMatPropConstCmd->SetToBeBroadcasted(false);

    fWorldMaterialCmd = new G4UIcmdWithAString("/opnovice2/detector/worldMaterial", this);
    fWorldMaterialCmd->SetGuidance("Set material of world.");
    fWorldMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fWorldMaterialCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fOpticalDir;
  delete fWorldMatPropVectorCmd;
  delete fWorldMatPropConstCmd;
  delete fWorldMaterialCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    //    FINISH
    if(command == fNPSAngleCmd)
    {
        fDetector->SetNPSAngle(fNPSAngleCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fNPSDistanceCmd)
    {
        fDetector->SetNPSDistance(fNPSDistanceCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fNPSShieldThicknessCmd)
    {
        fDetector->SetNPSShieldThickness(fNPSShieldThicknessCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fHCALAngleCmd)
    {
        fDetector->SetHCALAngle(fHCALAngleCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fHCALDistanceCmd)
    {
        fDetector->SetHCALDistance(fHCALDistanceCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fHCALShieldThicknessCmd)
    {
        fDetector->SetHCALShieldThickness(fHCALShieldThicknessCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fWindowThicknessCmd)
    {
        fDetector->SetWindowThickness(fWindowThicknessCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fTargetLengthCmd)
    {
        fDetector->SetTargetLength(fTargetLengthCmd->GetNewDoubleValue(newValue));
    }
    else if(command == fBeamlineCmd)
    {
        fDetector->SetBeamlineOn(fBeamlineCmd->GetNewIntValue(newValue));
    }
    else if(command == fUpdateCmd)
    {
        fDetector->UpdateGeometry();
    }
    else if(command == fWorldMatPropVectorCmd)
    {
        // Convert string to physics vector
        // string format is property name, then pairs of energy, value
        G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
        std::istringstream instring(newValue);
        G4String prop;
        instring >> prop;
        while(instring)
        {
            G4String tmp;
            instring >> tmp;
            if(tmp == "")
            {
                break;
            }
            G4double en = G4UIcommand::ConvertToDouble(tmp);
            instring >> tmp;
            G4double val = G4UIcommand::ConvertToDouble(tmp);
            mpv->InsertValues(en, val);
        }
        fDetector->AddWorldMPV(prop, mpv);
    }
    else if(command == fWorldMatPropConstCmd)
    {
        // Convert string to physics vector
        // string format is property name, then value
        // space delimited
        std::istringstream instring(newValue);
        G4String prop;
        G4String tmp;
        instring >> prop;
        instring >> tmp;
        G4double val = G4UIcommand::ConvertToDouble(tmp);
        fDetector->AddWorldMPC(prop, val);
    }
    else if(command == fWorldMaterialCmd)
    {
        fDetector->SetWorldMaterial(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
