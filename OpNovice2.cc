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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Initialized Modules
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

// Geant4 Libraries
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {
    // detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if(argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    // creates the Geant4 RunManager using the G4RunManagerFactory Factory class
    auto runManager = G4RunManagerFactory::CreateRunManager();

    DetectorConstruction* detCon = new DetectorConstruction(); // defines the Geometry and Detector of the Experiment
    runManager->SetUserInitialization(detCon);

    // Configures the Physics of the experiment
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    physicsList->RegisterPhysics(opticalPhysics);
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialization());

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // Get the Pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if(ui) // Interactive Mode
    {
        // Starting the Geant4 UI session for interactive mode
        UImanager->ApplyCommand("/control/execute ../macros/vis.mac");
        ui->SessionStart();
        delete ui; // Cleaning up the user interface session
    }
    else // Batch Mode
    {
        // Applying Geant4 control commands from the specified file
        G4String command  = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    // Cleaning up resources
    delete visManager; // Cleaning up the visual manager
    delete runManager; // Cleaning up the Geant4 RunManager
    return 0; // Exiting the program
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
