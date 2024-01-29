#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//---------------------------------------------------------------------------

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  :fPPhysicsList(pPhys)
{   
  fPhysDir = new G4UIdirectory("/opnovice2/physics/");
  fPhysDir->SetGuidance("physics control.");

  fPListCmd = new G4UIcmdWithAString("/opnovice2/physics/addPhysics",this);
  fPListCmd->SetGuidance("Add physics list (standard_opt3, QGSP_BIC_EMY, QGSP_BIS_HP, QGSP_BERT_HP, FTFP_BERT_HP");
}

//---------------------------------------------------------------------------

PhysicsListMessenger::~PhysicsListMessenger()
{;}

//---------------------------------------------------------------------------

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == fPListCmd )
   { fPPhysicsList->AddPhysicsList(newValue);}
}

//---------------------------------------------------------------------------

