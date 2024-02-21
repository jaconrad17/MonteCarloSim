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

// Standard C++ Libraries
#include <stdexcept>
#include <iostream>
#include <unistd.h>
#include <fstream>

// Geant4 Libraries
#include "globals.hh"
#include "CLHEP/Random/Random.h"
#include "G4RunManager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

// Initialized Modules
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "OutputManager.hh"
#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {
    // Setting Event Seed
    CLHEP::HepRandom::createInstance();
    
    // Generates a random seed for the event
    unsigned int seed = time(0) + (int)getpid();
    unsigned int devrandseed = 0;
    FILE *fdrand = fopen("/dev/urandom", "r");
    if (fdrand) {
        fread(&devrandseed, sizeof(int), 1, fdrand);
        seed += devrandseed;
        fclose(fdrand);
    }

    // Setting the seed for the random number generator
    CLHEP::HepRandom::setTheSeed(seed); 
    
    // detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if(argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }

    // creates the Geant4 RunManager using the G4RunManagerFactory Factory class
    auto runManager = G4RunManagerFactory::CreateRunManager();

    // creates a new instance of each class using dynamic memory
    PhysicsList* phys = new PhysicsList(); // defines the Physics of the Experiment
    DetectorConstruction* detCon = new DetectorConstruction(); // defines the Geometry and Detector of the Experiment
    PrimaryGeneratorAction* pga = new PrimaryGeneratorAction();
    RunAction* runAction = new RunAction(pga);
    OutputManager* outManager = new OutputManager(); // Outputs results to ROOT
    EventAction* event = new EventAction(runAction, outManager, pga);

    // Configuring the RunManager with user-defined classes
    runManager->SetUserInitialization(phys);
    runManager->SetUserInitialization(detCon);
    runManager->SetUserInitialization(pga);
    runManager->SetUserInitialization(runAction);
    runManager->SetUserAction(event);

    // Setting additional user actions
    runManager->SetUserAction(new SteppingAction(event));
    runManager->SetUserAction(new TrackingAction);

    // Initialize visualization
    G4VisManager* visManager = nullptr;

    // Get the Pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if(ui) // Interactive Mode
    {
        // Setting up Geant4 visualization manager for interactive mode
        visManager = new G4VisExecutive;
        visManager->Initialize();

        // Starting the Geant4 UI session for interactive mode
        ui->SessionStart();
        delete ui; // Cleaning up the user interface session
    }
    else // Batch Mode
    {
        if (argc > 1)
        {
            G4String fileName = argv[1];
            // Check if the specified file exists before applying the command
            if (std::ifstream(fileName))
            {
                // Applying Geant4 control commands from the specified file
                G4String command  = "/control/execute ";G4String command  = "/control/execute ";
                UImanager->ApplyCommand(command + fileName);

                // Running the simulation for the specified number of events
                G4String commandr = "/run/beamOn ";
                G4int nev = pga->GetNEvents();
                char snev[50];
                sprintf( snev, "%d", nev);
                UImanager->ApplyCommand(commandr+snev);
            }
            else
            {
                G4cerr << "Error: File" << fileName << " does not exist." << G4endl;
            }
        }
        else
        {
            G4cerr << "Error: No input file specified in batch mode." << G4endl;
        }
    }

    // Cleaning up resources
    delete outManager; // Cleaning up the output manager
    delete runManager; // Cleaning up the Geant4 RunManager
    return 0; // Exiting the program
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
