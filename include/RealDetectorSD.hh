#ifndef RealDetectorSD_h
#define RealDetectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4VHitsCollection.hh"
#include "RealDetectorHit.hh"

//---------------------------------------------------------------------------

class RealDetectorSD : public G4VSensitiveDetector
{
public:
    
    RealDetectorSD( G4String, G4int );
    ~RealDetectorSD();
    
    void   Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
    G4bool ProcessHits_constStep(const G4Step* astep,G4TouchableHistory* ROHist);
    void   EndOfEvent(G4HCofThisEvent*);
    
private:
    
    RealDetectorHitsCollection*  fCollection;
    G4int                    fDetID;
    G4int                    fNelements;
    G4int                    fNhits;
    G4int*                   fhitID;
    G4int*                   fHits;
    
};
#endif

//---------------------------------------------------------------------------






