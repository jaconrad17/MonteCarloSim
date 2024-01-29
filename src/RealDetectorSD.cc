#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4THitsCollection.hh"

#include "RealDetectorSD.hh"
#include "RealDetectorHit.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

RealDetectorSD::RealDetectorSD(G4String name, G4int)
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

RealDetectorSD::~RealDetectorSD()
{
}

//---------------------------------------------------------------------------

void RealDetectorSD::Initialize(G4HCofThisEvent* HCE)
{
    
    fCollection = new RealDetectorHitsCollection(SensitiveDetectorName,collectionName[0]);
    static G4int HCID = -1;
    
    if( HCID < 0 )  {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        
    }
    HCE->AddHitsCollection( HCID, (G4VHitsCollection*)fCollection );
}

//---------------------------------------------------------------------------

G4bool RealDetectorSD::ProcessHits( G4Step* aStep,G4TouchableHistory* )
{
    
    G4TouchableHistory*   theTouchable = (G4TouchableHistory*)(aStep->GetPostStepPoint()->GetTouchable());
    G4VPhysicalVolume*    volume       = theTouchable->GetVolume();
    G4int                 id           = volume->GetCopyNo();
    G4double              edep         = aStep->GetTotalEnergyDeposit();
    
    if (fhitID[id]==-1){
        
        RealDetectorHit* Hit = new RealDetectorHit;
        Hit->AddEnergy(edep);
        Hit->SetPrePosition(aStep->GetPreStepPoint()->GetPosition());
        Hit->SetPostPosition(aStep->GetPostStepPoint()->GetPosition());
        Hit->SetID(id);
        Hit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
        
        fhitID[id] = fCollection->insert(Hit) -1;
        fHits[fNhits++]=id;
    }
    else // This is not new
    {
        (*fCollection)[fhitID[id]]->AddEnergy(edep);
    }
    
    return true;
        
}

//---------------------------------------------------------------------------

G4bool RealDetectorSD::ProcessHits_constStep(const G4Step*,G4TouchableHistory*)
{
    return false;
}

//---------------------------------------------------------------------------

void RealDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
    
    for (G4int i=0;i<fNhits;i++) {
        fhitID[fHits[i]] = -1;
        fHits[i]         = 0;
    }
    fNhits = 0;
}

//---------------------------------------------------------------------------














