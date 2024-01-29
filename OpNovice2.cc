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
//
/// \file optical/OpNovice2/OpNovice2.cc
/// \brief Main program of the optical/OpNovice2 example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// from sFFG4MC
#include <stdexcept>
#include <iostream>
#include "globals.hh"
#include "CLHEP/Random/Random.h"

#include "G4RunManager.hh"
#include "G4UItcsh.hh"
#include "G4UIterminal.hh"

// from OpNovice2
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
    // Setting Event Seed
    CLHEP::HepRandom::createInstance();
    
    unsigned int seed = time(0) + (int) getpid();
    unsigned int devrandseed = 0;
    FILE *fdrand = fopen("/dev/urandom", "r");
    if ( fdrand ){
        fread(&devrandseed, sizeof(int), 1, fdrand);
        seed += devrandseed;
        fclose(fdrand);
    }
    
    CLHEP::HepRandom::setTheSeed(seed);
    
    // detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if(argc == 1) ui = new G4UIExecutive(argc, argv);

    auto runManager = G4RunManagerFactory::CreateRunManager();
    PhysicsList* phys = new PhysicsList();
    runManager->SetUserInitialization(phys);
    
    DetectorConstruction* detCon = new DetectorConstruction();
    runManager->SetUserInitialization(detCon);

    runManager->SetUserInitialization(new ActionInitialization());

    // initialize visualization
    G4VisManager* visManager = nullptr;
    /*G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();*/

    // get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if(ui)
    {
        // interactive mode
        //UImanager->ApplyCommand("/control/execute ../macros/vis.mac");
        visManager = new G4VisExecutive;
        visManager->Initialize();
        ui->SessionStart();
        delete ui;
    }
    else
    {
        // batch mode
        if (argc > 1)
        {
            G4String fileName = argv[1];
            // Check if the file exists before applying the command
            if (std::ifstream(fileName))
            {
                G4String command  = "/control/execute ";G4String command  = "/control/execute ";
                UImanager->ApplyCommand(command + fileName);
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

    // job termination
    delete visManager;
    delete outManager;
    delete runManager;
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
