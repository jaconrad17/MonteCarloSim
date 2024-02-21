#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4RunManager.hh"

/// Event action class
///

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction, Outputmanager* out, PrimaryGeneratorAction* pga);
    ~EventAction() override;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;
	
	void AddCount(G4int count) {fCounter += count; }
	void AddCount_left(G4int count_left) {fCounter_left += count_left; }
	void AddCount_right(G4int count_right) {fCounter_right += count_right; }
	void AddCount_mid(G4int count_mid) {fCounter_mid += count_mid; }

  private:
    RunAction* fRunAction = nullptr;
	G4int fCounter = 0;
	G4int fCounter_left = 0;
    G4int fCounter_right = 0;
	G4int fCounter_mid = 0;
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
