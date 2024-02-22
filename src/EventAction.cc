#include "EventAction.hh"
#include "RunAction.hh"
#include "VirtualDetectorHit.hh"
#include "RealDetectorHit.hh"
#include "OutputManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(OutputManager* out, PrimaryGeneratorAction* pga)
: fOutManager(out), fPGA(pga)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
    if( evt->GetEventID() == 0)
        fOutManager->InitOutput();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
    G4int event_id = evt->GetEventID();
    if( event_id%100 == 0)
    {
        printf("Event &8d/r", event_id);
        fflush(stdout);
    }
    
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    VirtualDetectorHit* hit;
    G4int nvirtualhits = 0;
    RealDetectorHit* rhit;
    G4int nrealhits = 0;
    G4int nprimhits = 0;
    
    if(HCE)
    {
        G4int CollSize = HCE->GetNumberOfCollections();
        G4int hci = 0;
        
        fOutManager->ZeroArray();
        
        for(G4int = 0; i < CollSize, i++)
        {
            VirtualDetectorHitsCollection* hc;
            while(!(hc = static_cast<VirtualDetectorHitsCollection*>(HCE->GetHC(hci++))));
            G4int hc_nhits = hc->entries();
            
            if(hc_nhits == 0) continue;
            
            // Fill output hit arrays for virtual detectors
            else if(hc->GetName().contains("VirtualDetector"))
            {
                for(G4int j = 0; j < hc_nhits; j++)
                {
                    hit = static_cast<VirtualDetectorHit*>( hc->GetHit(j));
                    fOutManager->SetVirtualPDef( (G4ParticleDefinition*) hit->GetPDef());
                    fOutManager->SetVirtualP3( (G4ThreeVector) hit->GetMom());
                    fOutManager->SetVirtualPosPre( (G4ThreeVector) hit->GetPosPre());
                    fOutManager->SetVirtualPosPost( (G4ThreeVector) hit->GetPosPost());
                    fOutManager->SetVirtualVertex( (G4ThreeVector) hit->GetVertex());
                    fOutManager->SetVirtualTime( (G4double) hit->GetTime());
                    fOutManager->SetVirtualEnergy( (G4double) hit->GetEdep());
                    fOutManager->SetVirtualID( (G4int) hit->GetID());
                    fOutManager->SetVirtualPID( (G4int) hit->GetParentID());
                    fOutManager->SetVirtualTID( (G4int) hit->GetTrackID());
                    
                    fOutManager->FillVirtualArray( nvirtualhits);
                    nvirtualhits++;
                }
            }
            
            // Fill output hit arrays for real detectors
            else if( hc->GetName().contains("RealDetector"))
            {
                for(G4int j = 0; j <hc_nhits; j++)
                {
                    rhit = static_cast<RealDetectorHit*>( hc->GetHit(j));
                    fOutManager->SetRealPosPre( (G4ThreeVector) rhit->GetPosPre());
                    fOutManager->SetRealPosPost( (G4ThreeVector) rhit->GetPosPost());
                    fOutManager->SetRealEdep( (G4double) rhit->GetEdep());
                    fOutManager->SetRealID( (G4int) rhit->GetID());
                    fOutManager->SetRealTime( (G4double) rhit->GetTime());
                    
                    fOutManager->FillRealArray(nrealhits);
                    nrealhits++;
                }
            }
            
        }
        
        // Fill output data from Primary GeneratorAction
        G4int nprim = fPGA->GetNPrimaryParticles();
        for(G4int j = 0; j < nprim; j++)
        {
            fOutManager->SetPrimaryDirection( (G4ThreeVector)fPGA->GetDirection(j));
            fOutManager->SetPrimaryEnergy( (G4double)fPGA->GetEnergy(j));
            fOutManager->SetPrimaryPDef( (G4ParticleDefinition*)fPGA->GetPrimPDef(j));
            fOutManager->SetPrimaryVertex( (G4ThreeVector)fPGA->GetVertex(j));
            
            fOutManager->FillPrimaryArray(nprimhits);
            nprimhits++;
        }
        
        // Only fill output root tree if there are hits in the event
        if(nvirtualhits != 0 || nrealhits != 0)
        {
            fOutManager->FillTree((G4double)fPGA->GetWeight(), (G4int)fPGA->GetFlag(), (G4int)fPGA->GetNEvents());
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
