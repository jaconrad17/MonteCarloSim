#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//---------------------------------------------------------------------------

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  : fPPhysicsList(pPhys)
{   
  // Creates a directory for physics commands
  fPhysDir = new G4UIdirectory("/opnovice2/physics/");
  fPhysDir->SetGuidance("physics control.");

  // Command to add physics list
  fPListCmd = new G4UIcmdWithAString("/opnovice2/physics/addPhysics", this);
  fPListCmd->SetGuidance("Add physics list (standard_opt4, standard_opt3, LowE_Livermore, LowE_Penelope, QGSP_BIC_EMY, QGSP_BIS_HP, QGSP_BERT_HP, FTFP_BERT_HP");
}

//---------------------------------------------------------------------------

PhysicsListMessenger::~PhysicsListMessenger()
{
  // Destructor
}

//---------------------------------------------------------------------------

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fPListCmd)
   { 
     // Attempt to add the specified physics list
     fPPhysicsList->AddPhysicsList(newValue);
   }
  else
  {
    // Handle unrecognized commands or provide feedback
    G4cerr << "PhysicsListMessenger::SetNewValue - Error: Unrecognized command." << G4endl;
  }
}

//---------------------------------------------------------------------------

