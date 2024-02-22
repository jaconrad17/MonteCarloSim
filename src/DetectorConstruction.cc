//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file optical/OpNovice2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(), fDetectorMessenger(nullptr)
{
    // Dimensions
    fExpHall_x = fExpHall_y = fExpHall_z = 10.75 *m;
    fNPSAngle                            = 15.5 *deg;
    fNPSDist                             = 301.0 *cm;
    fNPSShieldThick                      = 1.0 *cm;
    fHCALAngle                           = 42.5 *deg;
    fHCALDist                            = 427.0 *cm;
    fHCALShieldThick                     = 5.0 *cm;
    fSCWinThick                          = 0.050 *cm;
    fTarLength                           = 10.0 *cm;
    fBeamline                            = 0;

    // MaterialPropertiesTable Initialization
    fWorldMPT      = new G4MaterialPropertiesTable(); // World
    fPbWO4MPT      = new G4MaterialPropertiesTable(); // NPS Electron Arm
    fNPSshieldMPT  = new G4MaterialPropertiesTable(); // Electron Arm Shields
    fHCALscintMPT  = new G4MaterialPropertiesTable(); // HCAL Proton Arm
    fHCALeabsMPT   = new G4MaterialPropertiesTable(); // HCAL Proton Arm Shields
    fHodoscintMPT  = new G4MaterialPropertiesTable(); // Hadron Arm Hodoscope
    fHCALshieldMPT = new G4MaterialPropertiesTable(); // Hadron Arm Shield
    

    
  fTank_x = fTank_y = fTank_z = 1.0 * cm;

  fTank = nullptr;

  fTankMPT    = new G4MaterialPropertiesTable();
  fScintMPT = new G4MaterialPropertiesTable();
  
  fSurfaceMPT = new G4MaterialPropertiesTable();
  fSurfaceMPT2 = new G4MaterialPropertiesTable();

  fSurface2 = new G4OpticalSurface("Surface2");
  fSurface2->SetType(dielectric_dielectric);
  fSurface2->SetFinish(polished);
  fSurface2->SetModel(unified);

  fSurface = new G4OpticalSurface("Surface");
  fSurface->SetType(dielectric_dielectric);
  fSurface->SetFinish(polished);
  fSurface->SetModel(unified);

  const G4int NUM = 6;
  G4double pp[NUM] = {2.0*eV, 2.2*eV, 2.4*eV, 2.6*eV, 2.8*eV, 3.0*eV};
  G4double rindex[NUM] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
  G4double rindex2[NUM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  G4double rindex3[NUM] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double reflectivity[NUM] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3};

  G4double reflectivity2[NUM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  G4double tran2[NUM] = {0., 0., 0., 0., 0., 0.};

  G4double tran[NUM] = {0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
  G4double absorption[NUM] = {3.448*m, 4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m};
  fTankMPT->AddProperty("RINDEX", pp, rindex3, NUM);      // was rindex not rindex3
  //fTankMPT->AddProperty("ABSLENGTH", pp, absorption, NUM);
  //fSurfaceMPT->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
  //fSurfaceMPT->AddProperty("TRANSMITTANCE",pp,tran,NUM);

  fSurfaceMPT2->AddProperty("REFLECTIVITY", pp, reflectivity2, NUM);
  fSurfaceMPT2->AddProperty("TRANSMITTANCE", pp, tran2, NUM);

  fWorldMPT->AddProperty("RINDEX", pp, rindex2, NUM);
  fScintMPT->AddProperty("RINDEX", pp, rindex3, NUM);


  fSurface2->SetMaterialPropertiesTable(fSurfaceMPT2);

  fSurface->SetMaterialPropertiesTable(fSurfaceMPT);

  fTank_LV  = nullptr;
  fWorld_LV = nullptr;
  rect_mid_LV = nullptr;
  cone_LV = nullptr;
  rem_cyl_LV = nullptr;
  rem_cyl2_LV = nullptr;
  rem_cyl3_LV = nullptr;
  rem_cyl4_LV = nullptr;
  rec_box_LV = nullptr;

    // Initialize Logical Volumes
    fWorld_LV = nullptr; // World
    fPbWO4_LV = nullptr; // NPS Electron Arm
    fNPSshield_LV = nullptr; // Electron Arm Shields
    fHCALscint_LV = nullptr; // HCAL Proton Arm
    fHCALeabs_LV = nullptr; // HCAL Proton Arm Shield
    fHodoscint_LV = nullptr; // Hadron Arm Hodoscope
    fHCALshield_LV = nullptr; // Hadron Arm Shield

  bend_ang = 1*rad;
  bend_rad = 24.75*cm;
  

    
    
  fTankMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
  
  fScintMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");


    // Set Materials
    fWorldMaterial      = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"); // World
    fPbWO4Material      = G4NistManager::Instance()->FindOrBuildMaterial("G4_PbWO4"); // NPS Electron Arm
    fNPSshieldMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"); // Electron Arm Shields
    fHCALscintMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // HCAL Proton Arm
    fHCALeabsMaterial   = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"); // HCAL Proton Arm Shield
    fHodoscintMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // Hadron Arm Hodoscope
    fHCALshieldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"); // Hadron Arm Shield
    
    // Initialize DetectorMessenger
    fDetMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() 
{
    delete fPbWO4MPT;
    delete fNPSshieldMPT;
    delete fHCALscintMPT;
    delete fHCALeabsMPT;
    delete fHodoscintMPT;
    delete fHCALshieldMPT;
    delete fWorldMPT;
    delete fDetMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4int SDcount = 1;

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    fPbWO4_material->SetMaterialPropertiesTable(fPbWO4MPT);
    
    
  fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
  fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  fScintMaterial->SetMaterialPropertiesTable(fScintMPT);
  fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);

  G4double thick = 0.25*cm;
  G4double len = thick * 20;

  // ------------- Volumes --------------
  // The experimental Hall
  G4Box* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

  fWorld_LV = new G4LogicalVolume(world_box, fWorldMaterial, "World", 0, 0, 0);

  fWorld_PV = new G4PVPlacement(0, G4ThreeVector(), fWorld_LV, "World", 0, false, 0);

  G4RotationMatrix* Rot = new G4RotationMatrix();
  Rot->rotateY(90*deg);
  Rot->rotateX(90*deg - bend_ang);

  G4RotationMatrix* Rott = new G4RotationMatrix();
  Rott->rotateX(-bend_ang);

  G4RotationMatrix* Rotty = new G4RotationMatrix();
  Rotty->rotateY(-90.*deg);
  Rotty->rotateZ(90.*deg);

  G4RotationMatrix* Ro = new G4RotationMatrix();
  Ro->rotateX(90.*deg);

  //---------------------------------------------------------------------------
  // Create scattering chamber, target and exit beamline
  // Modified from original NPS simulation:
  // https://github.com/gboon18/HallC_NPS
  //---------------------------------------------------------------------------
  
  if(fBeamline)
      BuildBeamline();
  
  else
      BuildTarget();
    
    //--------------------------------------------------------------------------- 
    // Create "NPS" electron arm
    //--------------------------------------------------------------------------- 
  
    G4double fPbWO4_X = 20.*mm; 
    G4double fPbWO4_Y = 20.*mm; 
    G4double fPbWO4_Z = 200.*mm;

    G4Box* fPbWO4_solid = new G4Box("fPbWO4_solid", 0.5*fPbWO4_X, 0.5*fPbWO4_Y, 0.5*fPbWO4_Z);
    G4LogicalVolume* fPbWO4_LV = new G4LogicalVolume(fPbWO4_solid, fPbWO4Material, "fPbWO4_LV");

    G4double NPS_x, NPS_z, NPS_th, NPS_ph;
    G4double NPS_xprime, NPS_yprime, NPS_zprime;
    char stmp[50];

    for(int ix = 0; ix < fNPSNrow; ix++)
        {
            for(int iy = 0; iy < fNPSNcol; iy++) 
            {
                sprintf(stmp, "nps&d", SDcount);

                NPS_x = 0.0;
                NPS_y = -fPbWO4_Y + (ix*fPbWO4_Y);
                NPS_z = fNPSDist;

                NPS_th = fNPSAngle;

                NPS_ph = 0 + ((360./fNPSNcol * iy) *deg;

                NPS_yprime = NPS_y * std::cos(NPS_th) + NPS_z * std::sin(NPS_th);
                NPS_zprime = -NPS_y * std::sin(NPS_th) + NPS_z * std::cos(NPS_th);

                NPS_xprime = NPS_x * std::cos(NPS_ph) + NPS_yprime * std::sin(NPS_ph);
                NPS_yprime = NPS_x * std::sin(NPS_ph) + NPS_yprime * std::cos(NPS_ph);

                G4Transform3D fNPS_t3d = G4Translate3D(G4ThreeVector(NPS_xprime, NPS_yprime, NPS_zprime)) 
                  * G4RotateZ3D(NPS_ph).inverse() * G4RotateX3D(NPS_th).inverse();

                fDetVol[SDcount] new G4PVPlacement(fNPS_t3d, fPbWO4_LV, stmp, fWorld_LV, false, SDcount);

                SDcount++;
            }
        }

    //--------------------------------------------------------------------------- 
    // Create electron arm shields
    //---------------------------------------------------------------------------

    G4double fNPSshield_X = fPbWO4_X;
    G4double fNPSshield_Y = ((fNPSNrow) * fPbWO4_Y) + 5. *mm;
    G4double fNPSshield_Z = fNPSShieldThick;

    G4Box* fNPSshield_solid = new G4Box("NPSshield_solid", 0.5*fNPSshield_X, 0.5*fNPSshield_Y, 0.5*fNPSshield_Z);

    G4LogicalVolume* fNPSshield_LV = new G4LogicalVolume(fNPSshield_solid, fNPSshieldMaterial, "NPSshield_LV");

    G4double NPSshield_x, NPSshield_y, NPSshield_th, NPSshield_ph;
    G4double NPSshield_xprime, NPSshield_yprime, NPSshield_zprime;
    for(int iy = 0; iy < fNPSNcol; iy++)
    {
        sprintf(stmp, "npsshield%d", iy);

        NPSshield_x = 0.0;
        NPSshield_y = -fPbWO4_Y;
        NPSshield_z = fNPSDist - 0.5*fPbWO4_Z - fNPSshield_Z - 20 *mm;

        NPSshield_th = fNPSAngle;

        NPSshield_ph = 0 + ((360./fNPSNcol) * iy) *deg;

        NPSshield_yprime = NPSshield_z * std::sin(NPSshield_th); 
        NPSshield_zprime = NPSshield_z * std::cos(NPSshield_th); 
    
        NPSshield_xprime = NPSshield_x * std::cos(NPSshield_ph) + NPSshield_yprime * std::sin(NPSshield_ph); 
        NPSshield_yprime = NPSshield_x * std::sin(NPSshield_ph) + NPSshield_yprime * std::cos(NPSshield_ph); 

        G4Transform3D fNPSshield_t3d = G4Translate3D(G4ThreeVector(NPSshield_xprime, NPSshield_yprime, 0.0)) * G4RotateZ3D(NPSshield_ph).inverse() * G4Translate3D(G4ThreeVector(0.0, 0.0, NPSshield_zprime)) * G4RotateX3D(NPSshield_th).inverse();

        new G4PVPlacement(fNPSshield_t3d, fNPSshield_LV, stmp, fWorld_LV, false, 0);
    }

    //---------------------------------------------------------------------------
    // Create "HCAL" proton arm
    //--------------------------------------------------------------------------- 


  G4Tubs* rect_mid_curve = new G4Tubs("rect_mid1", bend_rad, bend_rad + (thick * 2), 2.49*cm, 0.*deg, bend_ang);
  G4Box* rect_mid_straight = new G4Box("rect_mid2", 2.49*fTank_x, thick, len);
  G4UnionSolid* rect_mid = new G4UnionSolid("rect_mid", rect_mid_curve, rect_mid_straight, Rot, G4ThreeVector((bend_rad + thick)*std::cos(bend_ang) - len*std::cos(90*deg - bend_ang), (bend_rad + thick)*std::sin(bend_ang) + len*std::sin(90*deg - bend_ang), 0));

  rect_mid_LV = new G4LogicalVolume(rect_mid, fTankMaterial, "rect_mid", 0, 0, 0);
  rect_mid_PV = new G4PVPlacement(Rotty, G4ThreeVector(0, (bend_rad + thick), 0), rect_mid_LV, "rect_mid", fWorld_LV, false, 0);


  // PMT
  G4Tubs* cone = new G4Tubs("Cone", 0., 2.75*cm, 5*cm, 0.*deg, 360.0*deg);
  cone_LV = new G4LogicalVolume(cone, fWorldMaterial, "Cone", 0, 0, 0);
  cone_PV = new G4PVPlacement(Rott, G4ThreeVector(0, (bend_rad + thick) * (1 - std::cos(bend_ang)) + 5*cm*std::sin(bend_ang) + 2*len*std::cos(90*deg - bend_ang), - (bend_rad + thick) * std::sin(bend_ang) - 5*cm*std::cos(bend_ang) - 2*len*std::sin(90*deg - bend_ang)), cone_LV, "Cone", fWorld_LV, false, 0);

  // scintillator
  G4Box* scint = new G4Box("Scint", 2.49*cm, thick, 5*cm);
  G4LogicalVolume* scint_LV = new G4LogicalVolume(scint, fTankMaterial, "Scint", 0, 0, 0);
  scint_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 5*cm), scint_LV, "Scint", fWorld_LV, false, 0);


    //---------------------------------------------------------------------------
    // Set Logical Attributes
    //---------------------------------------------------------------------------
    
    // Sensitive Detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    fVirtualDetectoerSD = new VirtualDetectorSD("VirtualDetectorSD", fNSD);
    SDman->AddNewDetector(fVirtualDetectorSD);
    pbwo4_log->SetSensitiveDetector(fVirtualDetectorSD);
    scintabs_log->SetSensitiveDetector(fVirtualDetectorSD);
    hodoscint_log->SetSensitiveDetector(fVirtualDetectorSD);

    fRealDetectorSD = new RealDetectorSD("RealDetectorSD", fNSD);
    SDman->AddNewDetector(fRealDetectorSD);
    pbwo4_log->SetSensitiveDetector(fRealDetectorSD);
    scintabs_log->SetSensitiveDetector(fRealDetectorSD);
    hodoscint_log->SetSensitiveDetector(fRealDetectorSD);
    
    // Visualization
    fLogicTarget->SetVisAttributes(G4Colour::Blue());
    pbwo4_LV->SetVisAttributes(G4Colour::Red());
    NPSshield_LV->SetVisAttributes(G4Colour::Cyan());
    scintabs_LV->SetVisAttributes(G4Colour::Yellow());
    hodoscint_LV->SetVisAttributes(G4Colour::Green());
    HCALshield_LV->SetVisAttributes(G4Colour::Cyan());

    feabs_log->SetVisAttributes(G4VisAttributes::Invisible);
    fWorld_LV->SetVisAttributes(G4VisAttributes::Invisible);
    f

    //---------------------------------------------------------------------------
    return fWorld_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfaceSigmaAlpha(G4double v)
{
  fSurface->SetSigmaAlpha(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfacePolish(G4double v)
{
  fSurface->SetPolish(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface polish set to: " << fSurface->GetPolish() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPV(const G4String& prop,
                                      G4MaterialPropertyVector* mpv)
{
  fTankMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const G4String& prop,
                                       G4MaterialPropertyVector* mpv)
{
  fWorldMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPV(const G4String& prop,
                                         G4MaterialPropertyVector* mpv)
{
  fSurfaceMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPC(const G4String& prop, G4double v)
{
  fTankMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPC(const G4String& prop, G4double v)
{
  fWorldMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPC(const G4String& prop, G4double v)
{
  fSurfaceMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fWorldMaterial != pmat)
  {
    fWorldMaterial = pmat;
    if(fWorld_LV)
    {
      fWorld_LV->SetMaterial(fWorldMaterial);
      fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "World material set to " << fWorldMaterial->GetName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTankMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fTankMaterial != pmat)
  {
    fTankMaterial = pmat;
    if(fTank_LV)
    {
      fTank_LV->SetMaterial(fTankMaterial);
      fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
      fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "Tank material set to " << fTankMaterial->GetName() << G4endl;
  }
}

void DetectorConstruction::SetBendRadius(G4double radius)
{
  bend_rad = radius*cm;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  world_PV = DetectorConstruction::Construct();
}

void DetectorConstruction::SetBendAngle(G4double angle)
{
  bend_ang = angle*deg;
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  world_PV = DetectorConstruction::Construct();
}
