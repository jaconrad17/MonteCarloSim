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
/// \file optical/OpNovice2/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

enum { EPGA_BEAM, EPGA_ROOT};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    
    virtual void GeneratePrimaries(G4Event*);
    void SetUpROOTInput (TString filename);
    
    void SetMode(G4int mo) {fMode = mo;}
    void SetNEvents(G4int nev) {fNevents = nev;}
    
    G4ParticleGun* GetParticleGun() { return fParticleGun; };
    
    G4int GetMode() {return fMode;}
    G4int GetNEvents() {return fNevents;}
    
    G4double GetWeight() {return (G4double)fWeight;}
    G4int GetFlag() {return fFlag;}
    G4int GetNPrimaryParticles() {return fNPrimParticles;}
    G4ThreeVector GetVertex(G4int i) {return G4ThreeVector(fVx[i], fVy[i], fVz[i]);}
    G4ThreeVector GetDirection(G4int i) {return G4ThreeVector(fPx[i], fPy[i], fPz[i]);}
    G4double GetEnergy(G4int i) {return fE[i];}
    G4ParticleDefinition* GetPrimPDef(G4int i){return fPDefinition[i];}
    
    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);
    void SetRandomDirection(G4bool val = true);
    G4bool GetPolarized() { return fPolarized; };
    G4double GetPolarization() { return fPolarization; }
    
private:
    PrimaryGeneratorMessenger* fGunMessenger;
    G4ParticleGun* fParticleGun;
    G4GeneralParticleSource* fParticleSource;
    G4ParticleTable* fParticleTable;
    
    G4int fMode;
    G4int fNevents;
    G4int fEvent;
    
    G4bool fRandomDirection;
    G4bool fPolarized;
    G4double fPolarization;
    
    static const G4int fMaxprim = 50;
    
    G4ParticleDefinitions* fPDefinition[fMaxprim];
    
    TFile* fGenFile;
    TTree* fGenTree;
    Int_t fNGenBranches;
    
    Int_t fNPrimParticles;
    Double_t fWeight;
    Int_t fFlag;
    
    Int_t fPDG[fMaxprim];
    
    Double_t fVx[fMaxprim];
    Double_t fVy[fMaxprim];
    Double_t fVz[fMaxprim];
    Double_t fPx[fMaxprim];
    Double_t fPy[fMaxprim];
    Double_t fPz[fMaxprim];
    Double_t fE[fMaxprim];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
