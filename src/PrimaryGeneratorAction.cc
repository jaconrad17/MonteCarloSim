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

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
  , fParticleGun(0)
{
    fMode = EPGA_BEAM;
    fNevents = -1;
    fEvent = 0;
    fWeight = 0;
    fFlag = 0;
    
    fGenFile = NULL;
    fGenTree = NULL;
    fNGenBranches = 0;
    
    G4int n_particle = 1;
    fParticleSource = new G4GeneralParticleSource();
    fParticleGun = new G4ParticleGun(n_particle);
    
    // create a messenger for this class
    fGunMessenger = new PrimaryGeneratorMessenger(this);
    
    fRandomDirection = false;
    fPolarized = false;
    fPolarization = 0.;
    // default kinematic
    //
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("e-");
    
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleTime(0.0 * ns);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0 * cm, 0.0 * cm, -69.0 * cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(2.5 * geV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    if(fGenFile){
        fGenFile->Close();
        delete fGenFile;
    }
    
    delete fParticleGun;
    delete fParticleSource;
    delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    switch (fMode)
    {
    case EPGA_BEAM:
        fParticleSource->GeneratePrimaryVertex(anEvent);

        fFlag = 999;
        fNPrimParticles = 1;
        fWeight = 60.e-6 / 1.602e-19;

        fVx[0] = fParticleSource->GetParticlePosition().getX();
        fVy[0] = fParticleSource->GetParticlePosition().getY();
        fVz[0] = fParticleSource->GetParticlePosition().getZ();
        fPx[0] = fParticleSource->GetParticleMomentumDirection().getX();
        fPy[0] = fParticleSource->GetParticleMomentumDirection().getY();
        fPz[0] = fParticleSource->GetParticleMomentumDirection().getZ();
        fE[0] = fParticleSource->GetParticleEnergy();
        fPDefinition[0] = fParticleSource->GetParticleDefinition();

        break;

    case EPGA_ROOT:
        if (fGenTree)
        {
            fGenTree->GetEvent(fEvent++);

            for (Int_t j = 0; j < fNPrimParticles; j++)
            {
                fPDefinition[j] = fParticleTable->FindParticle(fPDG[j]);

                fParticleGun->SetParticlePosition(G4ThreeVector(fVx[j] * cm, fVy[j] * cm, fVz[j] * cm));
                fParticleGun->SetParticleDefinition(fPDefinition[j]);
                fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fPx[j], fPy[j], fPz[j]).unit());
                fParticleGun->SetParticleEnergy(fE[j]);
                fParticleGun->GeneratePrimaryVertex(anEvent);

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

    default:
        G4cout << "Unknown mode given to PrimaryGeneratorAction (0 for BEAM, 1 for ROOT)" << G4endl;
    }
}

void PrimaryGeneratorAction::SetUpROOTInput(TString filename)
{
    fMode = EPGA_ROOT;
    
    fGenFile = new TFile(filename);
    if(!fGenFile)
        G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find filename" << G4endl;
    
    fGenTree = dynamic_cast<TTree*>(fGenFile->Get("TGen"));
    if(!fGenTree)
        G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find tree TGen" << G4endl;
    
    fNGenBranches = fGenTree->GetNbranches();
    TObjArray* objarray = fGenTree->GetListOfBranches();
    
    for(Int_t i = 0; i < fNGenBranches; i++)
    {
      
      TBranch *branch = dynamic_cast<TBranch*>     (objarray->At(i));
      TString  bname  = TString( const_cast<char*> (branch->GetName()));

      if(bname == "Nparticles") branch->SetAddress(&fNPrimParticles);
      if(bname == "weight") branch->SetAddress(&fWeight);
      if(bname == "flag") branch->SetAddress(&fFlag);

        if( fNPrimParticles < fMaxprim)
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

    G4int nev = fGenTree->GetEntries();
    if(fNevents == -1 || fNevents > nev)
      fNevents = nev;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  if(fParticleGun->GetParticleDefinition() !=
     G4OpticalPhoton::OpticalPhotonDefinition())
  {
    G4ExceptionDescription ed;
    ed << "The particleGun is not an opticalphoton.";
    G4Exception("PrimaryGeneratorAction::SetOptPhotonPolar", "OpNovice2_004",
                JustWarning, ed);
    return;
  }

  fPolarized    = true;
  fPolarization = angle;

  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2       = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if(modul2 > 0.)
    e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  G4ThreeVector polar =
    std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::SetRandomDirection(G4bool val)
{
  fRandomDirection = val;
}

