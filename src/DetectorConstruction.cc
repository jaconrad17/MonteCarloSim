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

    fPbWO4Material->SetMaterialPropertiesTable(fPbWO4MPT);
    fNPSshieldMaterial->SetMaterialPropertiesTable(fNPSshieldMPT);
    fHCALscintMaterial->SetMaterialPropertiesTable(fHCALscintMPT);
    fHCALeabsMaterial->SetMaterialPropertiesTable(fHCALeabsMPT);
    fHodoscintMaterial->SetMaterialPropertiesTable(fHodoscintMPT);
    fHCALshieldMaterial->SetMaterialPropertiesTable(fHCALshieldMPT);

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

                G4Transform3D fNPS_t3d = G4Translate3D(G4ThreeVector(NPS_xprime, NPS_yprime, NPS_zprime)) * G4RotateZ3D(NPS_ph).inverse() * G4RotateX3D(NPS_th).inverse();

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

    // 44 pairs of 10mm scint + 13mm Fe = 1012 mm
    G4int fHCALNpairs = 44;
    G4double fHCALeabs_Z = 13.0*mm;
    G4double fHCALscint_X = 150.0*mm;
    G4double fHCALscint_Y = 150.0*mm;
    G4double fHCALscint_Z = 1012.0*mm;

    G4Box* fHCALscint_solid = new G4Box("fHCALscint_solid", 0.5*fHCALscint_X, 0.5*fHCALscint_Y, 0.5*fHCALscint_Z);
    G4LogicalVolume* fHCALscint_LV = new G4LogicalVolume(fHCALscint_solid, fHCALscintMaterial, "fHCALscint_LV");
    
    G4Box* fHCALeabs_solid = new G4Box("fHCALeabs_solid", 0.5*fHCALscint_X, 0.5*fHCALscint_Y, 0.5*fHCALeabs_Z);
    G4LogicalVolume* fHCALeabs_LV = new G4LogicalVolume(fHCALeabs_solid, fHCALeabsMaterial, "fHCALeabs_LV");

    for(int iz = 0; iz < fHCALNpairs; iz++) 
    {
        sprintf(stmp, "feabs%d", iz);
        new G4PVPlacement(0, G4ThreeVector(0., 0., -fHCALscint_Z/2 + (iz+0.5)*(10*mm + fHCALeabs_Z)), fHCALeabs_LV, stmp, fHCALscint_LV, false, 9999);
    }

    G4double HCAL_x, HCAL_y, HCAL_z, HCAL_th, HCAL_ph;
    G4double HCAL_xprime, HCAL_yprime, HCAL_zprime;

    for(int ix = 0; ix < fHCALNrow; ix++)
    {
            for(int iy = 0; iy < fHCALNcol, iy++)
            {
                sprintf(stmp, "hcal%d", SDcount);
                HCAL_x = 0.0;
                HCAL_y = -fHCALscint_Y + (ix*fHCALscint_Y);
                HCAL_z = fHCALDist;

                HCAL_th = fHCALAngle;

                HCAL_ph = 0 + ((360./fHCALNcol) * iy) *deg;

                HCAL_yprime = HCAL_y * * std::cos(HCAL_th) + HCAL_z * std::sin(HCAL_th);
                HCAL_zprime = -HCAL_y * std::sin(HCAL_th) + HCAL_z * std::cos(HCAL_th); 

                HCAL_xprime = HCAL_x * std::cos(HCAL_ph) + HCAL_yprime * std::sin(HCAL_ph); 
                HCAL_yprime = HCAL_x * std::sin(HCAL_ph) + HCAL_yprime * std::cos(HCAL_ph); 

                G4Transform3D HCAL_t3d = G4Translate3D(G4ThreeVector(HCAL_xprime, HCAL_yprime, HCAL_zprime)) * G4RotateZ3D(HCAL_ph).inverse() * G4RotateX3D(HCAL_th).inverse();

                fDetVol[SDcount] = new G4PVPlacement(fHCAL_t3d, fHCALscint_LV, stmp, fWorld_LV, false, SDcount);

                SDcount++;
            }
    }

    //--------------------------------------------------------------------------- 
    // Create hadron arm hodoscope layer
    //--------------------------------------------------------------------------- 

    G4double fHodoscint_X = 30.0 *mm;
    G4double fHodoscint_Y = 30.0 *mm;
    G4double fHodoscint_Z = 100.0 *mm;

    G4Box* fHodoscint_solid = new G4Box("fHodoscint_solid", 0.5*fHodoscint_X, 0.5*fHodoscint_Y, 0.5*fHodoscint_Z);

    G4LogicalVolume* fHodoscint_LV = new G4LogicalVolume(fHodoscint_solid, fHodoscintMaterial, "fHodoscint_LV");

    G4double HODO_x, HODO_y, HODO_z, HODO_th, HODO_ph;
    G4double HODO_xprime, HODO_yprime, HODO_zprime;

    for(int ix = 0; ix < fHodoNrow; ix++)
        {
            for(int iy = 0; iy < fHodoNcol; iy++)
            {
                sprintf(stmp, "hodo%d", SDcount);

                HODO_x = 0.0;
                HODO_y = -(15*fHodoscint_Y)/2 + (ix*fHodoscint_Y);
                HODO_z = fHCALDist - 0.5*fHCALscint_Z - 0.5*fHodoscint_Z - 20 *mm;

                HODO_th = fHCALAngle;

                HODO_ph = 0 + ((360./fHodoNcol) * iy) *deg;

                HODO_yprime = HODO_y * std::cos(HODO_th) + HODO_z * std::sin(HODO_th); 
                HODO_zprime = -HODO_y * std::sin(HODO_th) + HODO_z * std::cos(HODO_th);

                HODO_xprime = HODO_x * std::cos(HODO_ph) + HODO_yprime * std::sin(HODO_ph);
                HODO_yprime = HODO_x * std::sin(HODO_ph) + HODO_yprime * std::cos(HODO_ph);

                G4Transform3D HODO_t3d = G4Translate3D(G4ThreeVector(HODO_xprime, HODO_yprime, HODO_zprime)) * G4RotateZ3D(HODO_ph).inverse() * G4RotateX3D(HODO_th).inverse();

                fDetVol[SDcount] = new G4PVPlacement(HODO_t3d, fHodoscint_LV, stmp, fWorld_LV, false, SDcount);

                SDcount++;
            }
        }

    //--------------------------------------------------------------------------- 
    // Create hadron arm shield
    //--------------------------------------------------------------------------- 

    G4double fHCALshield_X = fHCALscint_X;
    G4double fHCALshield_Y = ((fHCALNrow) * fHCALscint_Y);
    G4double fHCALshield_Z = fNPSShieldThick;

    G4Box* fHCALshield_solid = new G4Box("fHCALshield_solid", 0.5*fHCALshield_X, 0.5*fHCALshield_Y, 0.5*fHCALshield_Z);

    G4LogicalVolume* fHCALshield_LV = new G4LogicalVolume(fHCALshield_solid, fHCALshieldMaterial, "fHCALshield_LV");

    G4double HCALshield_x, HCALshield_y, HCALshield_z, HCALshield_th, HCALshield_ph;
    G4double HCALshield_xprime, HCALshield_yprime, HCALshield_zprime;

    for(int iy = 0; iy < fHCALNcol; iy++) 
    {
        sprintf(stmp, "hcalshield%d", iy);

        HCALshield_x = 0.0;
        HCALshield_y = -fHCALscint_Y;
        HCALshield_z = fHCALDist - 0.5*fHCALscint_Z - fHodoscint_Z - fHCALshield_Z - 20 *mm;

        HCALshield_th = fHCALAngle;

        HCALshield_ph = 0 + ((360./fHCALNcol) * iy) *deg;

        HCALshield_yprime = HCALshield_z * std::sin(HCALshield_th); 
        HCALshield_zprime = HCALshield_z * std::cos(HCALshield_th); 

        HCALshield_xprime = HCALshield_x * std::cos(HCALshield_ph) + HCALshield_yprime * std::sin(HCALshield_ph); 
        HCALshield_yprime = HCALshield_x * std::sin(HCALshield_ph) + HCALshield_yprime * std::cos(HCALshield_ph);

        G4Transform3D HCALshield_t3d = G4Translate3D(G4ThreeVector(HCALshield_xprime, HCALshield_yprime, HCALshield_zprime)) * G4RotateZ3D(HCALshield_ph).inverse() * G4RotateX3D(HCALshield_th).inverse();

        new G4PVPlacement(HCALshield_t3d, fHCALshield_LV, stmp, fWorld_LV, false, 0);
    }

    //---------------------------------------------------------------------------
    // Set Logical Attributes
    //---------------------------------------------------------------------------
    
    // Sensitive Detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    fVirtualDetectoerSD = new VirtualDetectorSD("VirtualDetectorSD", fNSD);
    SDman->AddNewDetector(fVirtualDetectorSD);
    fPbWO4_LV->SetSensitiveDetector(fVirtualDetectorSD);
    fHCALscint_LV->SetSensitiveDetector(fVirtualDetectorSD);
    fHodoscint_LV->SetSensitiveDetector(fVirtualDetectorSD);

    fRealDetectorSD = new RealDetectorSD("RealDetectorSD", fNSD);
    SDman->AddNewDetector(fRealDetectorSD);
    fPbWO4_LV->SetSensitiveDetector(fRealDetectorSD);
    fHCALscint_LV->SetSensitiveDetector(fRealDetectorSD);
    fHodoscint_LV->SetSensitiveDetector(fRealDetectorSD);
    
    // Visualization
    fLogicTarget->SetVisAttributes(G4Colour::Blue());
    fPbWO4_LV->SetVisAttributes(G4Colour::Red());
    fNPSshield_LV->SetVisAttributes(G4Colour::Cyan());
    fHCALscint_LV->SetVisAttributes(G4Colour::Yellow());
    fHodoscint_LV->SetVisAttributes(G4Colour::Green());
    fHCALshield_LV->SetVisAttributes(G4Colour::Cyan());

    fHCALeabs_LV->SetVisAttributes(G4VisAttributes::Invisible);
    fWorld_LV->SetVisAttributes(G4VisAttributes::Invisible);
    f

    //---------------------------------------------------------------------------
    return fWorld_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
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
void DetectorConstruction::AddWorldMPC(const G4String& prop, G4double v)
{
  fWorldMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
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
void DetectorConstruction::BuildBeamline()
{
    G4Material* VacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* ChamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* WindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
    G4Material* TargetCellMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* BeampipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    const G4double inch = 2.54*cm;

    G4double ChamberOuterRadius = 22.5*inch;
    G4double ChamberInnerRadius = 20.5*inch;

    G4double ChamberHeight = 24.25*2*inch;
    G4double WindowHeight = 19*inch;
    G4double WindowSubHeight = 15*inch;
    G4double WindowFrameHeight = 20*inch;
    G4double WindowClampHeight = 20*inch;

    G4double WindowFrameThickness = 1.25*inch;
    G4double WindowClampThickness = 0.750*inch;
    G4double WindowThickness = fSCWinThick;

    G4double TargetLength = fTarLength;
    G4double TargetRadius = 0.5*50.*mm;
    G4double TargetCellLength = (0.125 + 150.)*mm;
    G4double TargetCellRadius = (0.125 + 0.5*50.)*mm;
    G4double TargetWindowThickness = 0.125*mm;

    G4double WindowInnerJoint1OuterRadius = 0.5*1.469*inch;
    G4double WindowInnerJoint1Thickness = 0.5*inch;
    G4double WindowInnerJoint2InnerRadius = 0.5*1.068*inch;
    G4double WindowInnerJoint2OuterRadius = 0.5*1.50*inch;
    G4double WindowInnerJoint2Thickness  = 0.109*inch;

    G4double WindowOuterJoint2InnerRadius = 0.5*0.68*inch;
    G4double WindowOuterJoint2OuterRadius = 0.5*1.062*inch;
    G4double WindowOuterJoint2Thickness = (0.62 + 1.562 - 0.137 - 0.06)*inch;
    G4double WindowOuterJoint2_1InnerRadius = WindowOuterJoint2OuterRadius;
    G4double WindowOuterJoint2_1OuterRadius = 0.5*1.5*inch;
    G4double WindowOuterJoint2_1Thickness = 0.19*inch;
    G4double WindowOuterJoint2_1Position = 0.62*inch;

    G4double WindowOuterJoint2_2dx = 1.12*inch;
    G4double WindowOuterJoint2_2dy = 1.12*inch ;
    G4double WindowOuterJoint2_2dz = 0.137*inch;
    G4double WindowOuterJoint2_3InnerRadius = 0.5*0.68*inch;
    G4double WindowOuterJoint2_3OuterRadius = 0.5*0.738*inch;
    G4double WindowOuterJoint2_3Thickness = 0.06*inch;

    G4double Beampipe1Innerdx1 = 18.915*mm;
    G4double Beampipe1Innerdx2 = 63.4*mm;
    G4double Beampipe1Innerdy1 = Beampipe1Innerdx1;
    G4double Beampipe1Innerdy2 = Beampipe1Innerdx2;
    G4double Beampipe1Outerdx1 = 25.265*mm;
    G4double Beampipe1Outerdx2 = 69.749*mm;
    G4double Beampipe1Outerdy1 = Beampipe1Outerdx1;
    G4double Beampipe1Outerdy2 = Beampipe1Outerdx2;
    G4double Beampipe1Length = 1685.925*mm;
    G4double Beampipe2OuterRadius = 0.5*168.275*mm;
    G4double Beampipe2InnerRadius = 0.5*154.051*mm;
    G4double Beampipe2Length = 129.633*inch;

    G4double Beampipe2FrontCellThickness = Beampipe2OuterRadius -Beampipe2InnerRadius;
    G4double Beampipe3OuterRadius = 0.5*273.05*mm;
    G4double Beampipe3InnerRadius = 0.5*254.508*mm;
    G4double Beampipe3Length = 5102.225*mm;
    G4double Beampipe3FrontCellThickness = Beampipe3OuterRadius -Beampipe3InnerRadius;

    //---------------------------------------------------------------------------
    //Scattering Chamber
    //---------------------------------------------------------------------------

    G4RotationMatrix *xChambRot = new G4RotationMatrix;  
    xChambRot->rotateX(90*degree);                     
  
    G4Tubs* sChamberOuter = new G4Tubs("ChamberOuter_sol", ChamberInnerRadius, ChamberOuterRadius, 0.5*ChamberHeight, 0., twopi);

    G4double WindowStartTheta = (3.+10.+90.+180.)*pi/180.;
    G4double WindowDeltaTheta = 124.*pi/180.;
    G4double WindowFrameStartTheta = WindowStartTheta-(4.5)*pi/180.;
    G4double WindowFrameDeltaTheta = WindowDeltaTheta+9*pi/180.;

    G4Tubs* sWindowSub = new G4Tubs("WindowSub_sol", ChamberInnerRadius-1*cm, ChamberOuterRadius+1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);
  
    G4RotationMatrix *zWindowRot = new G4RotationMatrix;  
    zWindowRot->rotateZ(90*degree);                     
    G4SubtractionSolid* sChamber_sub_Window = new G4SubtractionSolid("Chamber_sub_Window", sChamberOuter, sWindowSub, zWindowRot, G4ThreeVector());

    G4LogicalVolume* LogicChamber = new G4LogicalVolume(sChamber_sub_Window, ChamberMaterial, "ChamberOuter_log");

    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicChamber, "ChamberOuter_pos", fWorld_LV, false, 0);

    //---------------------------------------------------------------------------
  
    G4Tubs* sWindowFrame_before_sub = new G4Tubs("WindowFrame_before_sub_sol", ChamberOuterRadius, ChamberOuterRadius + WindowFrameThickness, 0.5*WindowFrameHeight, WindowFrameStartTheta, WindowFrameDeltaTheta);

    G4Tubs* sWindowFrame_sub = new G4Tubs("WindowFrame_sub_sol", ChamberOuterRadius - 1*cm, ChamberOuterRadius + WindowFrameThickness + 1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);
  
    G4RotationMatrix *pseudoWindowRot = new G4RotationMatrix;  
    pseudoWindowRot->rotateZ(0*degree);                     
  
    G4SubtractionSolid* sWindowFrame = new G4SubtractionSolid("WindowFrame_sol", sWindowFrame_before_sub, sWindowFrame_sub, pseudoWindowRot, G4ThreeVector());

    G4LogicalVolume* LogicWindowFrame = new G4LogicalVolume(sWindowFrame, ChamberMaterial, "WindowFrame_log");

    G4RotationMatrix *zxWindowRot = new G4RotationMatrix(0,-90*degree,-90*degree);

    new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowFrame, "WindowFrame_pos", fWorld_LV, false, 0);

    //---------------------------------------------------------------------------

    G4Tubs* sWindowClamp_before_sub = new G4Tubs("WindowClamp_before_sub_sol", ChamberOuterRadius + WindowFrameThickness + WindowThickness, ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness, 0.5*WindowClampHeight, WindowFrameStartTheta, WindowFrameDeltaTheta);
  
    G4Tubs* sWindowClamp_sub = new G4Tubs("WindowClamp_sub_sol", ChamberOuterRadius + WindowFrameThickness + WindowThickness -1*cm, ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness +1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);

    G4SubtractionSolid* sWindowClamp = new G4SubtractionSolid("WindowClamp_sol", sWindowClamp_before_sub, sWindowClamp_sub, pseudoWindowRot, G4ThreeVector());
  
    G4LogicalVolume* LogicWindowClamp = new G4LogicalVolume(sWindowClamp, ChamberMaterial, "WindowClamp_log");

    new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowClamp, "WindowClamp_pos", fWorld_LV, false, 0);

    //---------------------------------------------------------------------------
    // Vacuum inside the chamber
    //---------------------------------------------------------------------------

    G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol", 0., ChamberInnerRadius, 0.5*ChamberHeight, 0., twopi);

    G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(sChamberInner, VacuumMaterial, "ChamberInner_log");

    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fWorld_LV, false, 0);

    //---------------------------------------------------------------------------
    // Target Cell & Target
    //---------------------------------------------------------------------------

    G4RotationMatrix *xTargetRot = new G4RotationMatrix;  
    xTargetRot->rotateX(-90*degree);                     

    G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol", 0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0.,twopi); 

    G4LogicalVolume* LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_log");   

    new G4PVPlacement(xTargetRot, G4ThreeVector(0., -43.9*cm, 0.), LogicTargetCell, "TargetCell_pos",  LogicInnerChamber, false, 0);

    //---------------------------------------------------------------------------

    G4Tubs* sTarget = new G4Tubs("Target_sol", 0., TargetRadius, 0.5*TargetLength, 0.,twopi); 

    fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_log");  
  
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0 ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::BuildTarget()
{
    G4Material* VacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Material* TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_lH2");
    G4Material* TargetCellMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    const G4double inch = 2.54*cm;

    G4double ChamberInnerRadius = 20.5*inch;
    G4double ChamberHeight = 3*24.25*2*inch;

    G4double BeamlineInnerRadius = 200. *mm;
    G4double BeamlineHeight = 1300. *cm;

    G4double TargetLength = fTarLength;
    G4double TargetRadius = 0.5*50.*mm;
    G4double TargetCellLength = (0.125 + 150.)*mm;
    G4double TargetCellRadius = (0.125 + 0.5*50.)*mm;
    G4double TargetWindowThickness = 0.125*mm;

    //---------------------------------------------------------------------------
    // Vacuum inside scattering chamber and exit beamline
    //---------------------------------------------------------------------------

    G4RotationMatrix *xChambRot = new G4RotationMatrix;
    xChambRot->rotateX(90*degree);

    G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol", 0., ChamberInnerRadius, 0.5*ChamberHeight, 0., twopi);
    G4Tubs* sBeamlineInner = new G4Tubs("BeamlineInner_sol", 0., BeamlineInnerRadius, BeamlineHeight, 0., twopi);
    
    G4UnionSolid* vacUnion = new G4UnionSolid("vacUnion", sChamberInner, sBeamlineInner, xChambRot, G4ThreeVector(0., -BeamlineHeight, 0.));
    G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(vacUnion, VacuumMaterial, "ChamberInner_LV")

    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fWorld_LV, false, 0);

    LogicInnerChamber->SetVisAttributes(G4VisAttributes::Invisible);

    //---------------------------------------------------------------------------
    // Target Cell & Target
    //---------------------------------------------------------------------------

    G4RotationMatrix *xTargetRot = new G4RotationMatrix;
    xTargetRot->rotateX(-90*degree);

    G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol", 0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0., twopi);

    G4LogicalVolume* LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_LV");

    new G4PVPlacement(xTargetRot, G4ThreeVector(), LogicTargetCell, "TargetCell_pos", LogicInnerChamber, false, 0);

    //---------------------------------------------------------------------------

    G4Tubs* sTarget = new G4Tubs("Target_sol", 0., TargetRadius, 0.5*TargetLength, 0., twopi);

    fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_LV");

    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0);
}
