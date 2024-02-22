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
/// \file optical/OpNovice2/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Include necessary header files
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4String.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "TObjArray.h"
#include "TBranch.h"
#include "TString.h"

// Use the CLHEP namespace
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Constructor for PrimaryGeneratorAction class
PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
  , fParticleGun(0)
{
  // Initialize member variables
  fMode = EPGA_BEAM; // Assume an beam setup
  fNevents = 1000; // Set a default number of events
  fEvent = 0;
  fWeight = 0;
  fFlag = 0;
    
  fGenFile = NULL;
  fGenTree = NULL;
  fNGenBranches = 0;

  // Set up a particle gun
  G4int n_particle = 1;
  fParticleSource = new G4GeneralParticleSource();
  fParticleGun = new G4ParticleGun(n_particle);
  fParticleTable = G4ParticleTable::GetParticleTable();

  // Set default values for particle gun
  fPolarized = false;
  fPolarization = 0.;

  // Create a messenger for this class to handle UI commands
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor for PrimaryGeneratorAction class
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  // Check if fGenFile is not NULL, close it and delete it
  if(fGenFile)
  {
    fGenFile->Close();
    delete fGenFile;
  }

  // Delete allocated memory for fParticleGun, fParticleSource, and fGunMessenger
  delete fParticleGun;
  delete fPartcleSource;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// GeneratePrimaries method for generating primary particles
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Switch statement based on the value of fMode  
  switch (fMode)
  {
    // Case for optical physics involving electrons
    case EPGA_BEAM: 
      
      if(fParticleGun->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
      {
        if(fPolarized)
          SetOptPhotonPolar(fPolarization);
        else
          SetOptPhotonPolar();
      }
      // Generate a primary vertex using fParticleSource
      fParticleSource->GeneratePrimaryVertex(anEvent);

      // Set values for specific variables
      fFlag = 999;
      fNPrimParticles = 1;
      fWeight = 60.e-6 / 1.602e-19;
        
      /* Retrieve properties of the generated particle and store them in arrays
      (position, momentum direction, enery, and particle definition)
      These arrays will be used later for ROOT mode */
      fVx[0] = fParticleSource->GetParticlePosition().getX();
      fVy[0] = fParticleSource->GetParticlePosition().getY();
      fVz[0] = fParticleSource->GetParticlePosition().getZ();
      fPx[0] = fParticleSource->GetParticleMomentumDirection().getX();
      fPy[0] = fParticleSource->GetParticleMomentumDirection().getY();
      fPz[0] = fParticleSource->GetParticleMomentumDirection().getZ();
      fE[0] = fParticleSource->GetParticleEnergy();
      fPDefinition[0] = fParticleSource->GetParticleDefinition();
    }
    break;
    // Case for ROOT mode
    case EPGA_ROOT:
      // If the generator tree (fGenTree) exists
      if (fGenTree)
      {
        // Get the event from the generator tree
        fGenTree->GetEvent(fEvent++);

        // Loop over the number of primary particles
        for (Int_t j = 0; j < fNPrimParticles; j++)
        {
          // Find the particle definition from the particle table using PDG code
          fPDefinition[j] = fParticleTable->FindParticle(fPDG[j]);

          // Set properties for the particle gun based on ROOT input
          fParticleGun->SetParticlePosition(G4ThreeVector(fVx[j] * cm, fVy[j] * cm, fVz[j] * cm));
          fParticleGun->SetParticleDefinition(fPDefinition[j]);
          fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fPx[j], fPy[j], fPz[j]).unit());
          fParticleGun->SetParticleEnergy(fE[j]);

          // Generate a primary vertex using the particle gun
          fParticleGun->GeneratePrimaryVertex(anEvent);

          // Update the arrays with the properties of the generated particle
          fVx[j] = fParticleGun->GetParticlePosition().getX();
          fVy[j] = fParticleGun->GetParticlePosition().getY();
          fVz[j] = fParticleGun->GetParticlePosition().getZ();
          fPx[j] = fParticleGun->GetParticleMomentumDirection().getX();
          fPy[j] = fParticleGun->GetParticleMomentumDirection().getY();
          fPz[j] = fParticleGun->GetParticleMomentumDirection().getZ();
          fE[j] = fParticleGun->GetParticleEnergy();
          fPDefinition[j] = fParticleGun->GetParticleDefinition();
        }
      }
      break;

    // Default case for an unknown mode
    default:
      G4cout << "Unknown mode given to PrimaryGeneratorAction (0 for BEAM, 1 for ROOT)" << G4endl;
  }
}

// Method to set up ROOT input for the generator
void PrimaryGeneratorAction::SetUpROOTInput(TString filename)
{
  // Set the mode to EPGA_ROOT
  fMode = EPGA_ROOT;

  // Open the ROOT file
  fGenFile = new TFile(filename);
  if(!fGenFile)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find filename" << G4endl;

  // Get the generator tree from the ROOT file
  fGenTree = dynamic_cast<TTree*>(fGenFile->Get("TGen"));
  if(!fGenTree)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find tree TGen" << G4endl;

  // Get the number of branches in the generator tree
  fNGenBranches = fGenTree->GetNbranches();
  TObjArray* objarray = fGenTree->GetListOfBranches();

  // Loop over the branches and set addresses for specific variable
  for(Int_t i = 0; i < fNGenBranches; i++)
  {
    TBranch *branch = dynamic_cast<TBranch*>     (objarray->At(i));
    TString  bname  = TString( const_cast<char*> (branch->GetName()));

    if(bname == "Nparticles") branch->SetAddress(&fNPrimParticles);
    if(bname == "weight") branch->SetAddress(&fWeight);
    if(bname == "flag") branch->SetAddress(&fFlag);

    if(fNPrimParticles < fMaxprim)
    {
      if(bname == "vx") branch->SetAddress(fVx);
      if(bname == "vy") branch->SetAddress(fVy);
      if(bname == "vz") branch->SetAddress(fVz);
      if(bname == "px") branch->SetAddress(fPx);
      if(bname == "py") branch->SetAddress(fPy);
      if(bname == "pz") branch->SetAddress(fPz);
      if(bname == "E") branch->SetAddress(fE);
      if(bname == "pdg") branch->SetAddress(fPDG);
    }
  }

  // Get the number of events in the generator tree
  G4int nev = fGenTree->GetEntries();
  if(fNevents == -1 || fNevents > nev)
    fNevents = nev;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Generate a random angle between 0 and 360 degrees
void PrimaryGeneratorAction::SetOptPhotonPolar()
{
  // Generate a random angle between 0 and 360 degrees
  G4double angle = G4UniformRand() * 360.0 * deg;
  // Call the SetOptPhotonPolar method with the generated angle
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// SetOptPhotonPolar method to set the polarization of optical photons wth a specified angle
void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  // Check if the particle gun is an optical photon
  if(fParticleGun->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
  {
    // If not an optical photon, issue a warning and return
    G4ExceptionDescription ed;
    ed << "The particleGun is not an opticalphoton.";
    G4Exception("PrimaryGeneratorAction::SetOptPhotonPolar", "OpNovice2_004",
                JustWarning, ed);
    return;
  }

  // Set polarization-related variables
  fPolarized    = true;
  fPolarization = angle;

  // Define vectors for calculations
  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2 = product * product;

  // Define perpendicular and parallel vectors based on polarization angle
  G4ThreeVector e_perpend(0., 0., 1.);
  if(modul2 > 0.)
    e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  // Calculate the polarization vector
  G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;

  // Set the particle polarization using the calculated vector
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

