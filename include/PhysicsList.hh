#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;

// Definition of the PhysicsList class derived from the G4VModularPhysicsList
class PhysicsList: public G4VModularPhysicsList
{
public:
// Constructor
PhysicsList();
// Destructor
virtual ~PhysicsList();

// Method to construct particles and their properties
void ConstructParticle();

// Method to set energy cuts for particles
void SetCuts();
void SetCutForGamma(G4double);
void SetCutForElectron(G4double);
void SetCutForPositron(G4double);

// Method to add a physics list by name
void AddPhysicsList(const G4String& name);

// Method to set the Verbose Level
void SetVerboseLevel(G4int level);

// Method to construct procsses related to physics
void ConstructProcess();

// Method to add a package (e.g., optical physics)
void AddPackage(const G4String& name);

// Method to construct photo-nuclear processes
void ConstructPhotoNuclear();
    
private:
// Configuration for electromagnetic physics
G4EmConfigurator em_config;

// Energy cuts for particles
G4double cutForGamma;
G4double cutForElectron;
G4double cutForPositron;

// Flags to check if certain physical processes are registered
G4bool helIsRegistered;
G4bool bicIsRegistered;
G4bool biciIsRegistered;
G4bool locIonIonInelasticIsRegistered;
G4bool radioactiveDecayIsRegistered;

// Names and instances of various physics constructors
G4String emName;
G4VPhysicsConstructor* emPhysicsList;
G4VPhysicsConstructor* decPhysicsList;
G4VPhysicsConstructor* raddecayList;

// Additional physics constructors for optical processes
G4VPhysicsConstructor* opticalPhysicsList;

// Vector to store instances of hadron physics constructors
std::vector<G4VPhysicsConstructor*> hadronPhys;

// Messenger class for handling UI commands related to the physics list
PhysicsListMessenger* pMessenger;    
};

#endif
