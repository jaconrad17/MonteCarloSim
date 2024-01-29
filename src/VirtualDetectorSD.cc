#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4THitsCollection.hh"

#include "VirtualDetectorSD.hh"
#include "VirtualDetectorHit.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

VirtualDetectorSD::VirtualDetectorSD(G4String name, G4int)
  :G4VSensitiveDetector(name)
{
    collectionName.insert(G4String("SDHits")+name);
    fCollection = NULL;
    fNelements  = 10000;
    fNhits      = 0;
    fhitID      = new G4int[fNelements];
    fHits       = new G4int[fNelements];
    for(G4int i=0; i<fNelements; i++) fhitID[i] = -1;
    for(G4int i=0; i<fNelements; i++) fHits[i]  = 0;
}

//---------------------------------------------------------------------------

VirtualDetectorSD::~VirtualDetectorSD()
{
}

//---------------------------------------------------------------------------

void VirtualDetectorSD::Initialize(G4HCofThisEvent* HCE)
{
    
    fCollection = new VirtualDetectorHitsCollection(SensitiveDetectorName,collectionName[0]);
    static G4int HCID = -1;
    
    if( HCID < 0 )  {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        
    }
    HCE->AddHitsCollection( HCID, (G4VHitsCollection*)fCollection );
}

//---------------------------------------------------------------------------

G4bool VirtualDetectorSD::ProcessHits( G4Step*,G4TouchableHistory* )
{
    return false;
}

//---------------------------------------------------------------------------

G4bool VirtualDetectorSD::ProcessHits_constStep(const G4Step* aStep,G4TouchableHistory*)
{
    
    G4Track*              aTrack       = aStep->GetTrack();
    G4String              ParticleName = aTrack->GetDefinition()->GetParticleName();
    G4double              ene          = aTrack->GetKineticEnergy();
    G4int                 pid          = aTrack->GetParentID();
    G4int                 tid          = aTrack->GetTrackID();
    G4TouchableHistory*   theTouchable = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
    G4VPhysicalVolume*    volume       = theTouchable->GetVolume();
    
    G4int  id;
    if( volume )
        id = volume->GetCopyNo();
    else {
        G4cout << "ERROR in virtualSD" << G4endl;
        id = 9999;
    }
    
    VirtualDetectorHit* Hit = new VirtualDetectorHit;
    Hit->SetEnergy(ene);
    Hit->SetMomentum(aTrack->GetMomentum());
    Hit->SetPrePosition(aStep->GetPreStepPoint()->GetPosition());
    Hit->SetPostPosition(aStep->GetPostStepPoint()->GetPosition());
    Hit->SetVertex(aTrack->GetVertexPosition());
    Hit->SetPDef(aTrack->GetDefinition());
    Hit->SetID(id);
    Hit->SetTrackID(tid);
    Hit->SetParentID(pid);
    Hit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
    fhitID[id] = fCollection->insert(Hit) -1;
    fHits[fNhits++]=id;
    return true;
        
}

//---------------------------------------------------------------------------

void VirtualDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
    
    for (G4int i=0;i<fNhits;i++) {
        fhitID[fHits[i]] = -1;
        fHits[i]         = 0;
    }
    fNhits = 0;    
}

//---------------------------------------------------------------------------














