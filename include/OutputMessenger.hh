#ifndef OutputMessenger_h
#define OutputMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------

class OutputManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//---------------------------------------------------------------------------

class OutputMessenger: public G4UImessenger
{
public:
    OutputMessenger(OutputManager* outMana);
    ~OutputMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    OutputManager* fOutputManager = nullptr;
    G4UIdirectory* fOutputDir = nullptr;
    G4UIcmdWithAString* fOutFileCmd = nullptr;
};
#endif

//---------------------------------------------------------------------------

