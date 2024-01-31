// Include the PhysicsList header and the PhysicsListMessenger header
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

// Include necessary Geant4 headers
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

// Include various Geant4 physics lists
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"
#include "G4PhotoNuclearProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4OpticalPhysics.hh"


//---------------------------------------------------------------------------

// Constructor implementation
PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    // Initialize Loss Table Manager
    G4LossTableManager::Instance();

    // Set default cut values for particles
    defaultCutValue = 0.1*mm;
    cutForGamma     = defaultCutValue;
    cutForElectron  = defaultCutValue;
    cutForPositron  = defaultCutValue;

    // Set flags for registering physics processes
    helIsRegistered  = false;
    bicIsRegistered  = false;
    biciIsRegistered = false;
    locIonIonInelasticIsRegistered = false;
    radioactiveDecayIsRegistered = false;

    // Set verbose level
    SetVerboseLevel(1);
    
    // EM physics
    emPhysicsList = new G4EmStandardPhysics_option4(1);
    emName = G4String("emstandard_opt4");
    
    //  Construct Optical Physics
    std::unique_ptr<G4OpticalPhysics> opticalPhysics = std::make_unique<G4OpticalPhysics>();
    emPhysicsList->RegisterPhysics(opticalPhysics.get());
    
    // Decay physics and all particles
    decPhysicsList = new G4DecayPhysics();
    raddecayList = new G4RadioactiveDecayPhysics();

    // Initialize PhysicsListMessenger
    pMessenger = new PhysicsListMessenger(this);
}

//---------------------------------------------------------------------------

// Destructor implementation
PhysicsList::~PhysicsList() = default;

//---------------------------------------------------------------------------

// Method to construct particles
void PhysicsList::ConstructParticle()
{
    decPhysicsList->ConstructParticle();
}

//---------------------------------------------------------------------------

// Method to construct physics processes
void PhysicsList::ConstructProcess()
{
    // Transportation
    AddTransportation();
    
    // Electromagnetic physics list
    emPhysicsList->ConstructProcess();
    em_config.AddModels();
    
    // Decay physics list
    decPhysicsList->ConstructProcess();
    raddecayList->ConstructProcess();
    
    // hadronic physics lists
    for(const auto& hadronPhysItem : hadronPhys) 
    {
        hadronPhysItem->ConstructProcess();
    }
}

//---------------------------------------------------------------------------

// Method to add a physics list by name
void PhysicsList::AddPhysicsList(const G4String& name)
{
    // Display a message indicating the physics list being added
    if (verboseLevel>1) 
    {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }

    // Check if the physics list to be added is the same as the current electromagnetic physics list
    if (name == emName) 
        return;

    // Check different physics list names and perform corresponding actions
    if (name == "standard_opt4") 
    {
        // Activate G4EmStandardPhysics_option4
        emName = name;
        delete emPhysicsList;
        emPhysicsList.reset(new G4EmStandardPhysics_option4());
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
    }
    else if (name == "standard_opt3") 
    {
        // Activate G4EmStandardPhysics_option3
        emName = name;
        emPhysicsList.reset(new G4EmStandardPhysics_option3());
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;
    } 
    else if (name == "LowE_Livermore") 
    {
        // Activate G4EmLivermorePhysics
        emName = name;
        emPhysicsList.reset(new G4EmLivermorePhysics());
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;
    } 
    else if (name == "LowE_Penelope") 
    {
        // Activate G4EmPenelopePhysics
        emName = name;
        emPhysicsList.reset(new G4EmPenelopePhysics());
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmPenelopePhysics" << G4endl; 
    }
    // Add QGSP_BIC physics and other relevant physics constructors
    else if (name == "QGSP_BIC_EMY" || name == "QGSP_BIC_HP" || name == "QGSP_BERT_HP" || name == "FTFP_BERT_HP") 
    {
        // Select QGSP_BIC_EMY, QGSP_BIC_HP, QGSP_BERT_HP, or FTFP_BERT_HP physics
        AddPhysicsList("emstandard_opt3");
        switch(name) 
        {
            case name == "QGSP_BIC_EMY":
                hadronPhys.push_back(std::make_unique<G4HadronPhysicsQGSP_BIC>());
                hadronPhys.push_back(std::make_unique<G4NeutronTrackingCut>());
                break;
            case name == "QGSP_BIC_HP": 
                hadronPhys.push_back(std::make_unique<G4HadronPhysicsQGSP_BIC_HP>()); 
                break;
            case name == "QGSP_BERT_HP": 
                hadronPhys.push_back(std::make_unique<G4HadronPhysicsQGSP_BERT_HP>()); 
                break;
            case name == "FTFP_BERT_HP": 
                hadronPhys.push_back(std::make_unique<G4HadronPhysicsFTFP_BERT_HP>()); 
                break;
            default: G4cout << "Error in QGSP_BIC & FTFP_BERT physics list condition" << G4endl;
        }
        hadronPhys.push_back(std::make_unique<G4EmExtraPhysics>());
        hadronPhys.push_back(std::make_unique<G4HadronElasticPhysics>());
        hadronPhys.push_back(std::make_unique<G4StoppingPhysics>());
        hadronPhys.push_back(std::make_unique<G4IonBinaryCascadePhysics>());
    } 
    else 
    {
        // Display an error message if the physics list is not recognized
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << " is not defined" << G4endl;
    }
}

//---------------------------------------------------------------------------

// Method to set particle cut values
void PhysicsList::SetCuts()
{
    // Display information about setting particle cuts if verboseLevel is greater than 0
    if (verboseLevel > 0)
    {
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
    }
    
    // Set cut values for gamma first, then for e-, and finally for e+
    // This ordering is important because some processes for e+/e- need cut values for gamma
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");
    
    // Set cuts for detector and display the cut values if verboseLevel is greater than 0
    if (verboseLevel>0) DumpCutValuesTable();
}

//---------------------------------------------------------------------------

// Method to set specific cut value for gamma particles
void PhysicsList::SetCutForGamma(G4double cut)
{
    cutForGamma = cut;
    // Set cut value for gamma particles using SetParticleCuts
    SetParticleCuts("gamma", cutForGamma);
}

//---------------------------------------------------------------------------

// Method to set specific cut value for electron particles
void PhysicsList::SetCutForElectron(G4double cut)
{
    cutForElectron = cut;
    // Set cut value for electron particles using SetParticleCuts
    SetParticleCuts("e-", cutForElectron);
}

//---------------------------------------------------------------------------

// Method to set specific cut value for positron particles
void PhysicsList::SetCutForPositron(G4double cut)
{
    cutForPositron = cut;
    // Set cut value for positron particles using SetParticleCuts
    SetParticleCuts("e+", cutForPositron);
}

//---------------------------------------------------------------------------

// Method to construct photo-nuclear process for gamma particles
void PhysicsList::ConstructPhotoNuclear()
{
    // Get the process manager for gamma particles
    G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();

    // Check if the process manager is available
    if (pManager)
    {
        // Create a new photo-nuclear process and a corresponding Bertini cascade interface
        G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
        G4CascadeInterface* bertini = new G4CascadeInterface();

        // Set the maximum energy for Bertini cascade to 10 GeV
        bertini->SetMaxEnergy(10*GeV);

        // Register the Bertini cascade interface with the photo-nuclear process
        process->RegisterMe(bertini);
        
        // Add the photo-nuclear process as a discrete process to the gamma particle's process manager
        pManager->AddDiscreteProcess(process);
    }
    else
    {
        G4cerr << "PhysicsList::ConstructPhotoNuclear - Error: Process manager for gamma particles not found." << G4endl;
    }
}

//---------------------------------------------------------------------------
