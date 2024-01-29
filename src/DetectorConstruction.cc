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

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

// New Libraries
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"

// From sFFG4MC
#include "VirtualDetectorSD.hh"
#include "RealDetectorSD.hh"

#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"

#include <fstream>

using namespace CLHEP;
using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
    fNistManager = G4NistManager::Instance();
    fDetMessenger = new DetectorMessenger(this);
    
    fExpHall_x = fExpHall_y = fExpHall_z = 10.75*m;
    fNPSAngle = 15.5*deg;
    fNPSDist = 301.0*cm;
    fNPSShieldThick = 1.0*cm;
    
    fHCALAngle = 42.5*deg;
    fHCALDist = 427.0*cm;
    fHCALShieldThick = 5.0*cm;
    
    fSCWinThick = 0.050*cm;
    fTarLength = 10.0*cm;
    fBeamline = 0;
    
    fWorldMPT = new G4MaterialPropertiesTable();
    fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() 
{
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
    
    // ------------ Surfaces -------------
    
    // Air
    G4MaterialPropertiesTable* matAirProperties = new G4MaterialPropertiesTable();
    matAirProperties->AddProperty("RINDEX", energy, rindex, numPoints);
    matAirProperties->AddProperty("ABSLENGTH", energy, absorption, numPoints);
    matAir->SetMaterialPropertiesTable(matAirProperties);
    
    // NPS Electron Arm
    G4MaterialPropertiesTable* matPbWO4Properties = new G4MaterialPropertiesTable();
    matPbWO4Properties->AddProperty("RINDEX", energy, rindex, numPoints);
    matPbWO4Properties->AddProperty("ABSLENGTH", energy, absorption, numPoints);
    matPbWO4->SetMaterialPropertiesTable(matPbWO4Properties);
    
    // HCAL Proton Arm
    G4MaterialPropertiesTable* matFeProperties = new G4MaterialPropertiesTable();
    matFeProperties->AddProperty("RINDEX", energy, rindex, numPoints);
    matFeProperties->AddProperty("ABSLENGTH", energy, absorption, numPoints);
    matFe->SetMaterialPropertiesTable(matFeProperties);
    
    G4MaterialPropertiesTable* matScintillatorProperties = new G4MaterialPropertiesTable();
    matScintillatorProperties->AddConstProperty("SCINTILLATIONYIELD", 5000./MeV);
    matScintillatorProperties->AddConstProperty("RESOLUTIONSCALE", 1.0);
    matScintillatorProperties->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
    matScintillatorProperties->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
    matScintillatorProperties->AddConstProperty("YIELDRATIO", 1.0);
    matScintillator->SetMaterialPropertiesTable(matScintillatorProperties);
    
    // Hadron Arm Hodoscope Layer
    G4MaterialPropertiesTable* matHodoscintProperties = new G4MaterialPropertiesTable();
    matHodoscintProperties->AddConstProperties("SCINTILLATIONYIELD", 10000./MeV);
    matHodoscintProperties->AddConstProperties("RESOLUTIONSCALE", 1.0);
    matHodoscintProperties->AddConstProperties("FASTTIMECONSTANT", 1.*ns);
    matHodoscintProperties->AddConstProperties("SLOWTIMECONSTANT", 10.*ns);
    matHodoscintProperties->AddConstProperties("YIELDRATIO", 1.0);
    matHodoscint->SetMaterialPropertiesTable(matHodoscintProperties);
    
    // Hadron Arm Shield
    G4MaterialPropertiesTable* matHCALshieldProperties = new G4MaterialPropertiesTable();
    matHCALshieldProperties->AddProperty("RINDEX", energy, rindex, numPoints);
    matHCALshieldProperties->AddProperty("ABSLENGTH", energy, absorption, numPoints);
    matHCALshield->SetMaterialPropertiesTable(matHCALshieldProperties);

    // ------------- Volumes --------------
    // The experimental Hall
    G4Box* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

    fWorld_LV = new G4LogicalVolume(world_box, fWorldMaterial, "World", 0, 0, 0);

    G4VPhysicalVolume* world_PV = new G4PVPlacement(0, G4ThreeVector(), fWorld_LV, "World", 0, false, 0);

    //---------------------------------------------------------------------------
    // Create scattering chamber, target and exit beamline
    // Modified from original NPS simulation:
    // https://github.com/gboon18/HallC_NPS
    //---------------------------------------------------------------------------
    
    if( fBeamline)
        BuildOpticalBeamline();
    else
        BuildOpticalTarget();
    
    //---------------------------------------------------------------------------
    // Create "NPS" electron arm
    //---------------------------------------------------------------------------

    G4double pbwo4_X = 20.*mm;
    G4double pbwo4_Y = 20.*mm;
    G4double pbwo4_Z = 200.*mm;
    
    G4Box* pbwo4_solid = new G4Box("pbwo4_solid", 0.5*pbwo4_X, 0.5*pbwo4_Y, 0.5*pbwo4_Z);
    
    G4LogicalVolume* pbwo4_log = new G4LogicalVolume(pbwo4_solid,
                             fNistManager->FindOrBuildMaterial("G4_PbWO4"),
                             "pbwo4_log");
    
    G4double NPS_x, NPS_y, NPS_z, NPS_th, NPS_ph;
    G4double NPS_xprime, NPS_yprime, NPS_zprime;
    char stmp[50];
    
    for( int ix=0 ; ix < fNPSNrow ; ix++ ) {
        for( int iy=0 ; iy < fNPSNcol ; iy++ ) {
            
            sprintf( stmp, "nps%d", SDcount );
            
            NPS_x  = 0.0;
            NPS_y  = -pbwo4_Y + (ix*pbwo4_Y);
            NPS_z  = fNPSDist;
            
            NPS_th = fNPSAngle;
            
            NPS_ph = 0 + ((360./fNPSNcol) * iy) *deg;
            
            NPS_yprime = NPS_y * std::cos(NPS_th) + NPS_z * std::sin(NPS_th);
            NPS_zprime = -NPS_y * std::sin(NPS_th) + NPS_z * std::cos(NPS_th);
            
            NPS_xprime = NPS_x * std::cos(NPS_ph) + NPS_yprime * std::sin(NPS_ph);
            NPS_yprime = NPS_x * std::sin(NPS_ph) + NPS_yprime * std::cos(NPS_ph);
            
            G4Transform3D NPS_t3d = G4Translate3D(G4ThreeVector(NPS_xprime, NPS_yprime, NPS_zprime))
            * G4RotateZ3D(NPS_ph).inverse() * G4RotateX3D(NPS_th).inverse();
            
            fDetVol[SDcount] = new G4PVPlacement(NPS_t3d, pbwo4_log, stmp, fWorld_LV, false, SDcount);
            
            pbwo4_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
            
            G4OpticalSurface* NPSopticalSurface = new G4OpticalSurface("OpticalSurface");
            NPSopticalSurface->SetType(dielectric_dielectric);
            NPSopticalSurface->SetFinish(polished);
            NPSopticalSurface->SetModel(unified);
            
            new G4LogicalSkinSurface("Surface", pbwo4_log, NPSopticalSurface);
            
            G4Scintillation* NPSscintillationProcess = new G4Scintillation();
            pbwo4_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
            
            G4OpAbsorption* NPSabsorptionProcess = new G4OpAbsorption();
            pbwo4_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
            
            G4OpBoundaryProcess* NPSboundaryProcess = new G4OpBoundaryProcess();
            pbwo4_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
            
            SDcount++;
        }
    }
    
    //---------------------------------------------------------------------------
    // Create electron arm shields
    //---------------------------------------------------------------------------

    G4double NPSshield_X = pbwo4_X;
    G4double NPSshield_Y = ((fNPSNrow) * pbwo4_Y) + 5.*mm;
    G4double NPSshield_Z = fNPSShieldThick;
    
    G4Box* NPSshield_solid = new G4Box("NPSshield_solid", 0.5*NPSshield_X, 0.5*NPSshield_Y, 0.5*NPSshield_Z);
    
    G4LogicalVolume* NPSshield_log = new G4LogicalVolume( NPSshield_solid,
                                  fNistManager->FindOrBuildMaterial("G4_Pb"),
                                  "NPSshield_log");

    G4double NPSshield_x, NPSshield_y, NPSshield_z, NPSshield_th, NPSshield_ph;
    G4double NPSshield_xprime, NPSshield_yprime, NPSshield_zprime;
    
    for( int iy=0 ; iy < fNPSNcol ; iy++ ) {
      
      sprintf( stmp, "npsshield%d", iy );
      
      NPSshield_x  = 0.0;
      NPSshield_y  = -pbwo4_Y;
      NPSshield_z  = fNPSDist - 0.5*pbwo4_Z - NPSshield_Z - 20 *mm;
      
      NPSshield_th = fNPSAngle;
      
      NPSshield_ph = 0 + ((360./fNPSNcol) * iy) *deg;
      
      NPSshield_yprime = NPSshield_z * std::sin(NPSshield_th);
      NPSshield_zprime = NPSshield_z * std::cos(NPSshield_th);
      
      NPSshield_xprime = NPSshield_x * std::cos(NPSshield_ph) + NPSshield_yprime * std::sin(NPSshield_ph);
      NPSshield_yprime = NPSshield_x * std::sin(NPSshield_ph) + NPSshield_yprime * std::cos(NPSshield_ph);
      
      G4Transform3D NPSshield_t3d = G4Translate3D(G4ThreeVector(NPSshield_xprime, NPSshield_yprime, 0.0)) * G4RotateZ3D(NPSshield_ph).inverse()
        * G4Translate3D(G4ThreeVector(0.0, 0.0, NPSshield_zprime)) * G4RotateX3D(NPSshield_th).inverse();
      
      new G4PVPlacement(NPSshield_t3d, NPSshield_log, stmp, fWorld_LV, false, 0);
        
        NPSshield_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        G4OpticalSurface* NPSShieldOpticalSurface = new G4OpticalSurface("OpticalSurface");
        NPSShieldOpticalSurface->SetType(dielectric_metal);
        NPSShieldOpticalSurface->SetFinish(polished);
        NPSShieldOpticalSurface->SetModel(unified);
        
        new G4LogicalSkinSurface("Surface", NPSshield_log, NPSShieldOpticalSurface);
        
        G4OpAbsorption* NPSShieldabsorptionProcess = new G4OpAbsorption();
        NPSshield_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        G4OpBoundaryProcess* NPSShieldboundaryProcess = new G4OpBoundaryProcess();
        NPSshield_log->SetUserLimits(new G4UserLImits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
      
    }
    
    //---------------------------------------------------------------------------
    // Create "HCAL" proton arm
    //---------------------------------------------------------------------------

    // 44 pairs of 10mm scint + 13mm Fe = 1012 mm
    G4int    HCALNpairs = 44;
    G4double feabs_Z    = 13.0*mm;
    G4double scintabs_X = 150.0*mm;
    G4double scintabs_Y = 150.0*mm;
    G4double scintabs_Z = 1012.0*mm;
    
    G4Box* scintabs_solid = new G4Box("scintabs_solid", 0.5*scintabs_X, 0.5*scintabs_Y, 0.5*scintabs_Z);
    
    G4Box* feabs_solid    = new G4Box("feabs_solid", 0.5*scintabs_X, 0.5*scintabs_Y, 0.5*feabs_Z);
    
    G4LogicalVolume* scintabs_log = new G4LogicalVolume(scintabs_solid,
                                fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"),
                                "scintabs_log");
    
    G4LogicalVolume* feabs_log = new G4LogicalVolume(feabs_solid,
                             fNistManager->FindOrBuildMaterial("G4_Fe"),
                             "feabs_log");
    
    for( int iz=0 ; iz < HCALNpairs ; iz++ ) {
      
      sprintf( stmp, "feabs%d", iz );
      new G4PVPlacement(0, G4ThreeVector( 0., 0., -scintabs_Z/2 + (iz+0.5)*(10*mm + feabs_Z)  ),
                feabs_log, stmp, scintabs_log, false, 9999 );
    }
    
    G4double HCAL_x, HCAL_y, HCAL_z, HCAL_th, HCAL_ph;
    G4double HCAL_xprime, HCAL_yprime, HCAL_zprime;
    
    for( int ix=0 ; ix < fHCALNrow ; ix++ ) {
      for( int iy=0 ; iy < fHCALNcol ; iy++ ) {
        
        sprintf(stmp, "hcal%d", SDcount);
        
        HCAL_x  = 0.0;
        HCAL_y  = -scintabs_Y + (ix*scintabs_Y);
        HCAL_z  = fHCALDist;
        
        HCAL_th = fHCALAngle;
        
        HCAL_ph = 0 + ((360./fHCALNcol) * iy) *deg;
        
        HCAL_yprime = HCAL_y * std::cos(HCAL_th) + HCAL_z * std::sin(HCAL_th);
        HCAL_zprime = -HCAL_y * std::sin(HCAL_th) + HCAL_z * std::cos(HCAL_th);
        
        HCAL_xprime = HCAL_x * std::cos(HCAL_ph) + HCAL_yprime * std::sin(HCAL_ph);
        HCAL_yprime = HCAL_x * std::sin(HCAL_ph) + HCAL_yprime * std::cos(HCAL_ph);
        
        G4Transform3D HCAL_t3d = G4Translate3D(G4ThreeVector(HCAL_xprime, HCAL_yprime, HCAL_zprime))
      * G4RotateZ3D(HCAL_ph).inverse() * G4RotateX3D(HCAL_th).inverse();
        
        fDetVol[SDcount] = new G4PVPlacement(HCAL_t3d, scintabs_log, stmp, fWorld_LV, false, SDcount);
          
          scintabs_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpticalSurface* HCALopticalSurface = new G4OpticalSurface("OpticalSurface");
          HCALopticalSurface->SetType(dielectric_dielectric);
          HCALopticalSurface->SetFinish(polished);
          HCALopticalSurface->SetModel(unified);
          
          new G4LogicalSkinSurface("Surface", scintabs_log, HCALopticalSurface);
          
          G4Scintillation* HCALscintillationProcess = newG4Scinillation();
          scintabs_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpAbsorption* HCALabsorptionProcess = new G4OpAbsorption();
          scintabs_logs->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpBoundaryProcess* HCALboundaryProcess = new G4OpBoundaryProcess();
          scintabs_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        SDcount++;
      }
    }
    
    //---------------------------------------------------------------------------
    // Create hadron arm hodoscope layer
    //---------------------------------------------------------------------------

    G4double hodoscint_X = 30.0 *mm;
    G4double hodoscint_Y = 30.0 *mm;
    G4double hodoscint_Z = 100.0*mm;
    
    G4Box* hodoscint_solid = new G4Box("hodoscint_solid", 0.5*hodoscint_X, 0.5*hodoscint_Y, 0.5*hodoscint_Z);
    
    G4LogicalVolume* hodoscint_log = new G4LogicalVolume(hodoscint_solid,
                                 fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"),
                                 "hodoscint_log");
    
    G4double HODO_x, HODO_y, HODO_z, HODO_th, HODO_ph;
    G4double HODO_xprime, HODO_yprime, HODO_zprime;

    for( int ix=0 ; ix < fHodoNrow ; ix++ ) {
      for( int iy=0 ; iy < fHodoNcol ; iy++ ) {
        
        sprintf( stmp, "hodo%d", SDcount );
        
        HODO_x  = 0.0;
        HODO_y  = -(15*hodoscint_Y)/2 + (ix*hodoscint_Y);
        HODO_z  = fHCALDist - 0.5*scintabs_Z - 0.5*hodoscint_Z - 20 *mm;
        
        HODO_th = fHCALAngle;
        
        HODO_ph = 0 + ((360./fHodoNcol) * iy) *deg;
        
        HODO_yprime = HODO_y * std::cos(HODO_th) + HODO_z * std::sin(HODO_th);
        HODO_zprime = -HODO_y * std::sin(HODO_th) + HODO_z * std::cos(HODO_th);
        
        HODO_xprime = HODO_x * std::cos(HODO_ph) + HODO_yprime * std::sin(HODO_ph);
        HODO_yprime = HODO_x * std::sin(HODO_ph) + HODO_yprime * std::cos(HODO_ph);
        
        G4Transform3D HODO_t3d = G4Translate3D(G4ThreeVector(HODO_xprime, HODO_yprime, HODO_zprime))
      * G4RotateZ3D(HODO_ph).inverse() * G4RotateX3D(HODO_th).inverse();
        
        fDetVol[SDcount] = new G4PVPlacement(HODO_t3d, hodoscint_log, stmp, fWorld_LV, false, SDcount);
          
          hodoscint_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpticalSurface* HODOopticalSurface = new G4OpticalSurface("OpticalSurface");
          HODOopticalSurface->SetType(dielectric_dielectric);
          HODOopticalSurface->SetFinish(polished);
          HODOopticalSurface->SetModel(unified);
          
          new G4LogicalSkinSurface("Surface", hodoscint_log, HODOopticalSurface);
          
          G4Scintillation* HODOscintillationProcess = new G4Scintillation();
          hodoscint_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpAbsorption* HODOabsorptionProcess = new G4OpAbsorption();
          hodoscint_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
          
          G4OpBoundaryProcess* HODOboundaryProcess = new G4OpBoundaryProcess();
          hodoscint_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        SDcount++;
      }
    }
    
    //---------------------------------------------------------------------------
    // Create hadron arm shield
    //---------------------------------------------------------------------------

    G4double HCALshield_X = scintabs_X;
    G4double HCALshield_Y = ((fHCALNrow) * scintabs_Y) ;
    G4double HCALshield_Z = fNPSShieldThick;
    
    G4Box* HCALshield_solid = new G4Box("HCALshield_solid", 0.5*HCALshield_X, 0.5*HCALshield_Y, 0.5*HCALshield_Z);
    
    G4LogicalVolume* HCALshield_log = new G4LogicalVolume(HCALshield_solid,
                              fNistManager->FindOrBuildMaterial("G4_Pb"),
                              "HCALshield_log");
    
    G4double HCALshield_x, HCALshield_y, HCALshield_z, HCALshield_th, HCALshield_ph;
    G4double HCALshield_xprime, HCALshield_yprime, HCALshield_zprime;
    
    for( int iy=0 ; iy < fHCALNcol ; iy++ ) {
      
      sprintf( stmp, "hcalshield%d", iy );
      
      HCALshield_x  = 0.0;
      HCALshield_y  = -scintabs_Y;
      HCALshield_z  = fHCALDist - 0.5*scintabs_Z - hodoscint_Z - HCALshield_Z - 20 *mm;
      
      HCALshield_th = fHCALAngle;
      
      HCALshield_ph = 0 + ((360./fHCALNcol) * iy) *deg;
      
      HCALshield_yprime = HCALshield_z * std::sin(HCALshield_th);
      HCALshield_zprime = HCALshield_z * std::cos(HCALshield_th);
      
      HCALshield_xprime = HCALshield_x * std::cos(HCALshield_ph) + HCALshield_yprime * std::sin(HCALshield_ph);
      HCALshield_yprime = HCALshield_x * std::sin(HCALshield_ph) + HCALshield_yprime * std::cos(HCALshield_ph);
      
      G4Transform3D HCALshield_t3d = G4Translate3D(G4ThreeVector(HCALshield_xprime, HCALshield_yprime, HCALshield_zprime))
        * G4RotateZ3D(HCALshield_ph).inverse() * G4RotateX3D(HCALshield_th).inverse();
      
      new G4PVPlacement(HCALshield_t3d, HCALshield_log, stmp, fWorld_LV, false, 0);
        
        HCALshield_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        G4OpticalSurface* HCALshieldOpticalSurface = new G4OpticalSurface("OpticalSurface");
        HCALshieldOpticalSurface->SetType(dielectric_metal);
        HCALshieldOpticalSurface->SetFinish(polished);
        HCALshieldOpticalSurface->SetModel(glisur);
        
        new G4LogicalSkinSurface("Surface", HCALshield_log, HCALshieldOpticalSurface);
        
        G4OpAbsorption* HCALshieldAbsorptionProcess = new G4OpAbsorption();
        HCALshield_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
        
        G4OpBoundaryProcess* HCALshieldBoundaryProcess = new G4OpBoundaryProcess();
        HCALshield_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX, 0., 0.));
      
    }
      
    //---------------------------------------------------------------------------
    // Set Logical Attributes
    //---------------------------------------------------------------------------

    // Senstive detector
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    fVirtualDetectorSD = new VirtualDetectorSD("VirtualDetectorSD", fNSD );
    SDman->AddNewDetector( fVirtualDetectorSD );
    pbwo4_log->SetSensitiveDetector( fVirtualDetectorSD );
    scintabs_log->SetSensitiveDetector( fVirtualDetectorSD );
    hodoscint_log->SetSensitiveDetector( fVirtualDetectorSD );

    fRealDetectorSD = new RealDetectorSD("RealDetectorSD", fNSD );
    SDman->AddNewDetector(fRealDetectorSD);
    pbwo4_log->SetSensitiveDetector(fRealDetectorSD);
    scintabs_log->SetSensitiveDetector(fRealDetectorSD);
    hodoscint_log->SetSensitiveDetector(fRealDetectorSD);

    // Visualisation
    fLogicTarget->SetVisAttributes(G4Colour::Blue());
    pbwo4_log->SetVisAttributes(G4Colour::Red());
    NPSshield_log->SetVisAttributes(G4Colour::Cyan());
    scintabs_log->SetVisAttributes(G4Colour::Yellow());
    hodoscint_log->SetVisAttributes(G4Colour::Green());
    HCALshield_log->SetVisAttributes(G4Colour::Cyan());

    feabs_log->SetVisAttributes(G4VisAttributes::Invisible);
    fWorld_LV->SetVisAttributes(G4VisAttributes::Invisible);

    //---------------------------------------------------------------------------

  return world_PV;
}

void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void DetectorConstruction::BuildOpticalBeamline()
{
    G4Material* VacuumMaterial = fNistManager->FindOrBuildMaterial("G4_Galactic");
    G4Material* ChamberMaterial = fNistManager->FindOrBuildMaterial("G4_Al");
    G4Material* WindowMaterial = fNistManager->FindOrBuildMaterial("G4_Al");
    G4Material* TargetMaterial = fNistManager->FindOrBuildMaterial("G4_lH2");
    G4Material* TargetCellMaterial = fNistManager->FindOrBuildMaterial("G4_Al");
    G4Material* BeampipeMaterial = fNistManager->FindOrBuildMaterial("G4_Al");
    
    // Setting Optical Properties for Aluminum Materials
    G4MaterialPropertiesTable* chamberMPT = new G4MaterialPropertiesTable();
    G4MaterialPropertiesTable* windowMPT = new G4MaterialPropertiesTable();
    G4MaterialPropertiesTable* targetCellMPT = new G4MaterialPropertiesTable();
    G4MaterialPropertiesTable* beampipeMPT = new G4MaterialPropertiesTable();
    
    const G4int numEntriesAl = 2;
    G4double photonEnergyAl[numEntriesAl]= {1.0*eV, 10.0*eV}; // Need to set
    G4double refractiveIndexAl[numEntriesAl] = {1.5, 1.5}; // Need to set
    
    G4double electronEnergy = 6.6 * GeV;
    G4double electronConductivityAl[1] = (1.0e-15 * siemens / cm);
    G4double electronConductivityH2[1] = (1.0e-15 * siemens / cm);
    
    chamberMPT->AddProperty("ELECTRICCONDUCTIVITY", electronEnergy, electronConductivityAl, 1);
    windowMPT->AddProperty(("ELECTRICCONDUCTIVITY", electronEnergy, electronConductivityAl, 1);
    targetCellMPT->AddProperty("ELECTRICCONDUCTIVITY", electronEnergy, electronConductivityAl, 1);
    beampipeMPT->AddProperty("ELECTRICCONDUCTIVITY", electronEnergy, electronConductivityAl, 1);
    
    // Setting Optical Properties for Liquid Hydrogen
    G4MaterialPropertiesTable* targetMPT = new G4MaterialPropertiesTable();
    
    const G4int numEntriesH2 = 2;
    G4double photonEnergyH2[numEntriesH2] = {1.0*eV, 10.0*eV};
    G4double refractiveIndexH2[numEntriesH2] = {1.3, 1.3};
    
    targetMPT->AddProperty("ELECTRICCONDUCTIVITY", electronEnergy, electronConductivityH2, 1);
    
    // Assign Optical Properties to Materials
    ChamberMaterial->SetMaterialPropertiesTable(chamberMPT);
    WindowMaterial->SetMaterialPropertiesTable(windowMPT);
    TargetMaterial->SetMaterialPropertiesTable(targetMPT);
    TargetCellMaterial->SetMaterialPropertiesTable(targetCellMPT);
    BeampipeMaterial->SetMaterialPropertiesTable(beampipeMPT);
    
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
    G4double WindowThickness = fSCWinThick; //0.020*inch;
    
    G4double TargetLength = fTarLength //100.*mm //(129.759-27.759+50.)*mm;
    G4double TargetRadius = 0.5*50.*mm //20.179*mm;
    G4double TargetCellLength = (0.125+150.)*mm; //(5.*129.759-29.759+50.)*mm;
    G4double TargetCellRadius = (0.125+0.5*50.)*mm; //(20.179+5.)*mm; //20.32*mm;
    G4double TargetWindowThickness = 0.125*mm; //5.*mm; //0.128*mm;
    
    G4double WindowInnerJoint1OuterRadius   = 0.5*1.469*inch;
    G4double WindowInnerJoint1Thickness     = 0.5*inch;
    G4double WindowInnerJoint2InnerRadius   = 0.5*1.068*inch;
    G4double WindowInnerJoint2OuterRadius   = 0.5*1.50*inch;
    G4double WindowInnerJoint2Thickness     = 0.109*inch;
    G4double WindowOuterJoint2InnerRadius   = 0.5*0.68*inch;
    G4double WindowOuterJoint2OuterRadius   = 0.5*1.062*inch;
    G4double WindowOuterJoint2Thickness     = (0.62 + 1.562 - 0.137 - 0.06)*inch;
    G4double WindowOuterJoint2_1InnerRadius = WindowOuterJoint2OuterRadius;
    G4double WindowOuterJoint2_1OuterRadius = 0.5*1.5*inch;
    G4double WindowOuterJoint2_1Thickness   = 0.19*inch;
    G4double WindowOuterJoint2_1Position    = 0.62*inch;
    G4double WindowOuterJoint2_2dx          = 1.12*inch;
    G4double WindowOuterJoint2_2dy          = 1.12*inch ;
    G4double WindowOuterJoint2_2dz          = 0.137*inch;
    G4double WindowOuterJoint2_3InnerRadius = 0.5*0.68*inch;
    G4double WindowOuterJoint2_3OuterRadius = 0.5*0.738*inch;
    G4double WindowOuterJoint2_3Thickness   = 0.06*inch;
    
    G4double Beampipe1Innerdx1    = 18.915*mm;
    G4double Beampipe1Innerdx2    = 63.4*mm;
    G4double Beampipe1Innerdy1    = Beampipe1Innerdx1;
    G4double Beampipe1Innerdy2    = Beampipe1Innerdx2;
    G4double Beampipe1Outerdx1    = 25.265*mm;
    G4double Beampipe1Outerdx2    = 69.749*mm;
    G4double Beampipe1Outerdy1    = Beampipe1Outerdx1;
    G4double Beampipe1Outerdy2    = Beampipe1Outerdx2;
    G4double Beampipe1Length      = 1685.925*mm;
    G4double Beampipe2OuterRadius = 0.5*168.275*mm;
    G4double Beampipe2InnerRadius = 0.5*154.051*mm;
    G4double Beampipe2Length      = 129.633*inch;//20180417 changed from 2620.138*mm;<-this is the actual length of the pipe.
    //however, there is a gap between beampipe1 and 2. In order to fill the gap. The length of the beampipe2 currently is longer.
    G4double Beampipe2FrontCellThickness = Beampipe2OuterRadius - Beampipe2InnerRadius; //made up
    G4double Beampipe3OuterRadius = 0.5*273.05*mm;
    G4double Beampipe3InnerRadius = 0.5*254.508*mm;
    G4double Beampipe3Length      = 5102.225*mm;
    G4double Beampipe3FrontCellThickness = Beampipe3OuterRadius - Beampipe3InnerRadius; //made up
    
    //---------------------------------------------------------------------------
    //Scattering Chamber
    //---------------------------------------------------------------------------

    G4RotationMatrix *xChambRot = new G4RotationMatrix;
    xChambRot->rotateX(90*degree);
    
    G4Tubs* sChamberOuter = new G4Tubs("ChamberOuter_sol", ChamberInnerRadius, ChamberOuterRadius, 0.5*ChamberHeight, 0., twopi);
    
    G4double WindowStartTheta = (3.+10.+90.+180.)*pi/180.;//window beginning position
    G4double WindowDeltaTheta = 124.*pi/180.;//window size 127deg
    G4double WindowFrameStartTheta = WindowStartTheta-(4.5)*pi/180.;//to match the center of the window
    G4double WindowFrameDeltaTheta = WindowDeltaTheta+9*pi/180.;//total 133 deg.

    G4Tubs* sWindowSub = new G4Tubs("WindowSub_sol", ChamberInnerRadius-1*cm, ChamberOuterRadius+1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);
    
    G4RotationMatrix *zWindowRot = new G4RotationMatrix;
    zWindowRot->rotateZ(90*degree);
    G4SubtractionSolid* sChamber_sub_Window = new G4SubtractionSolid("Chamber_sub_Window", sChamberOuter, sWindowSub, zWindowRot, G4ThreeVector());

    G4LogicalVolume* LogicChamber = new G4LogicalVolume(sChamber_sub_Window, ChamberMaterial, "ChamberOuter_log");
    LogicChamber->SetMaterialPropertiesTable(chamberMPT);

    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicChamber, "ChamberOuter_pos", fWorld_LV, false, 0);
    
    G4OpticalSurface* WindowSurface = new G4OpticalSurface("WindowSurface");
    WindowSurface->SetType(dielectric_dielectric);
    WindowSurface->SetFinish(polished);
    WindowSurface->SetModel(unified);
    
    new G4LogicalSkinSurface("WindowSurface", LogicChamber, WindowSurface);
    WindowSurface->SetMaterialPropertiesTable(windowMPT);
    
    //---------------------------------------------------------------------------

    G4Tubs* sWindowFrame_before_sub = new G4Tubs("WindowFrame_before_sub_sol",
                             ChamberOuterRadius, ChamberOuterRadius + WindowFrameThickness, 0.5*WindowFrameHeight,
                             WindowFrameStartTheta, WindowFrameDeltaTheta);

    G4Tubs* sWindowFrame_sub = new G4Tubs("WindowFrame_sub_sol", ChamberOuterRadius - 1*cm, ChamberOuterRadius + WindowFrameThickness + 1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);
    
    G4RotationMatrix *pseudoWindowRot = new G4RotationMatrix;
    pseudoWindowRot->rotateZ(0*degree);
    
    G4SubtractionSolid* sWindowFrame = new G4SubtractionSolid("WindowFrame_sol", sWindowFrame_before_sub, sWindowFrame_sub, pseudoWindowRot, G4ThreeVector());

    G4LogicalVolume* LogicWindowFrame = new G4LogicalVolume(sWindowFrame, ChamberMaterial, "WindowFrame_log");
    LogicWindowFrame->SetMaterialPropertiesTable(windowMPT);

    G4RotationMatrix *zxWindowRot = new G4RotationMatrix(0,-90*degree,-90*degree);

    new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowFrame, "WindowFrame_pos", fWorld_LV, false, 0);
    
    //---------------------------------------------------------------------------

    G4Tubs* sWindowClamp_before_sub = new G4Tubs("WindowClamp_before_sub_sol", ChamberOuterRadius + WindowFrameThickness + WindowThickness, ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness, 0.5*WindowClampHeight, WindowFrameStartTheta, WindowFrameDeltaTheta);
    
    G4Tubs* sWindowClamp_sub = new G4Tubs("WindowClamp_sub_sol", ChamberOuterRadius + WindowFrameThickness + WindowThickness -1*cm, ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness +1*cm, 0.5*WindowSubHeight, WindowStartTheta, WindowDeltaTheta);

    G4SubtractionSolid* sWindowClamp = new G4SubtractionSolid("WindowClamp_sol", sWindowClamp_before_sub, sWindowClamp_sub, pseudoWindowRot, G4ThreeVector());
    
    G4LogicalVolume* LogicWindowClamp = new G4LogicalVolume(sWindowClamp, ChamberMaterial, "WindowClamp_log");
    LogicWindowClamp->SetMaterialPropertiesTable(windowMPT);

    new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowClamp, "WindowClamp_pos", fWorld_LV, false, 0 );
    
    //---------------------------------------------------------------------------
    // Vacuum inside the chamber
    //---------------------------------------------------------------------------

    G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol", 0., ChamberInnerRadius, 0.5*ChamberHeight, 0., twopi);

    G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(sChamberInner, VacuumMaterial, "ChamberInner_log");

    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fWorld_LV, false, 0 );
    
    //---------------------------------------------------------------------------
    // Target Cell & Target
    //---------------------------------------------------------------------------

    G4RotationMatrix *xTargetRot = new G4RotationMatrix;
    xTargetRot->rotateX(-90*degree);

    G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol", 0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0.,twopi);

    G4LogicalVolume* LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_log");
    LogicTargetCell->SetMaterialPropertiesTable(targetCellMPT);

    new G4PVPlacement(xTargetRot, G4ThreeVector(0., -43.9*cm, 0.), LogicTargetCell, "TargetCell_pos",  LogicInnerChamber, false, 0 );
    
    //---------------------------------------------------------------------------

    G4Tubs* sTarget = new G4Tubs("Target_sol", 0., TargetRadius, 0.5*TargetLength, 0.,twopi);

    fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_log");
    fLogicTarget->SetMaterialPropertiesTable(targetMPT);
    
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0 );
    
}
void DetectorConstruction::BuildOpticalTarget()
{
    G4Material* VacuumMaterial = fNistManager->FindOrBuildManager("G4_Galactic");
    G4Material* TargetMaterial = fNistManager->FindOrBuildManager("G4_lH2");
    G4Material* TargetCellMaterial = fNistManager->FindOrBuildMaterial("G4_Al");
    
    // Setting Optical Properties for Aluminum Materials
    G4MaterialPropertiesTable* targetCellMPT = new G4MaterialPropertiesTable();
    
    const G4int numEntriesAl = 2;
    G4double photonEnergyAl[numEntriesAl]= {1.0*eV, 10.0*eV}; // Need to set
    G4double refractiveIndexAl[numEntriesAl] = {1.5, 1.5}; // Need to set
    
    targetCellMPT->AddProperty("RINDEX", photonEnergyAl, refractiveIndexAl, numEntriesAl);
    
    // Setting Optical Properties for Liquid Hydrogen
    G4MaterialPropertiesTable* targetMPT = new G4MaterialPropertiesTable();
    
    const G4int numEntriesH2 = 2;
    G4double photonEnergyH2[numEntriesH2] = {1.0*eV, 10.0*eV};
    G4double refractiveIndexH2[numEntriesH2] = {1.3, 1.3};
    
    targetMPT->AddProperty("RINDEX", photonEnergyH2, refractiveIndexH2, numEntriesH2);
    
    // Assign Optical Properties to Materials
    TargetMaterial->SetMaterialPropertiesTable(targetMPT);
    TargetCellMaterial->SetMaterialPropertiesTable(targetCellMPT);
    
    const G4double inch = 2.54*cm;
    
    G4double ChamberInnerRadius = 20.5*inch;
    G4double ChamberHeight = 3*24.25*2*inch;
    
    G4double BeamlineInnerRadius = 200.*mm;
    G4double BeamlineHeight = 1300.*cm;
    
    G4double TargetLength = fTarLength;
    G4double TargetRadius = 0.5*50.*mm;
    G4double TargetCellLength = (0.125+150.)*mm;
    G4double TargetWindowThickness = 0.125*mm;
    
    //---------------------------------------------------------------------------
    // Vacuum inside scattering chamber and exit beamline
    //---------------------------------------------------------------------------

    G4RotationMatrix *xChambRot = new G4RotationMatrix;
    xChambRot->rotateX(90*degree);

    G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol", 0., ChamberInnerRadius, 0.5*ChamberHeight, 0., twopi);

    G4Tubs* sBeamlineInner = new G4Tubs("BeamlineInner_sol", 0., BeamlineInnerRadius, BeamlineHight, 0., twopi);

    G4UnionSolid* vacUnion = new G4UnionSolid("vacUnion", sChamberInner, sBeamlineInner, xChambRot, G4ThreeVector(0., -BeamlineHight, 0.));

    G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(vacUnion, VacuumMaterial, "ChamberInner_log");
    
    new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fexpHall_log, false, 0 );

    LogicInnerChamber->SetVisAttributes(G4VisAttributes::Invisible);

    //---------------------------------------------------------------------------
    // Target Cell & Target
    //---------------------------------------------------------------------------

    G4RotationMatrix *xTargetRot = new G4RotationMatrix;
    xTargetRot->rotateX(-90*degree);

    G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol", 0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0.,twopi);

    G4LogicalVolume*   LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_log");

    new G4PVPlacement(xTargetRot, G4ThreeVector(), LogicTargetCell, "TargetCell_pos",  LogicInnerChamber, false, 0 );

    //---------------------------------------------------------------------------

    G4Tubs* sTarget = new G4Tubs("Target_sol",
                     0., TargetRadius, 0.5*TargetLength, 0.,twopi);

    fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_log");
    
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0 );
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
