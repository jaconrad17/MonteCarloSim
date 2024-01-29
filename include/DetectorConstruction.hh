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
/// \file optical/OpNovice2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4OpticalSurface.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "DetectorMessenger.hh"

//class DetectorMessenger;
class G4VPhysicalVolume;
class VirtualDetectorSD;
class RealDetectorSD;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    G4VPhysicalVolume* Construct();
    
    void UpdateGeometry();
    void BuildOpticalBeamLine();
    void BuildOpticalTarget();
    
    inline G4VPhysicalVolume* GetWorldPV() { return world_PV; };
    inline G4VPhysicalVolume* GetDetVol(G4int i) { return fDetVol[i]; };
    inline VirtualDetectorSD* GetVirtualDetectorSD() { return fVirtualDetectorSD; };
    inline G4int GetNoSD() { return fNSD; };
    
    void SetNPSAngle(G4double a) {fNPSAngle = a;  }
    void SetNPSDistance(G4double d) {fNPSDist = d;  }
    void SetNPSShieldThickness(G4double nt) {fNPSShieldThick = nt; }
    void SetHCALAngle(G4double an) {fHCALAngle = an; }
    void SetHCALDistance(G4double d) {fHCALDist = di; }
    void SetHCALShieldThickness(G4double t) {fHCALShieldThick = t;  }
    void SetWindowThickness(G4double th) {fSCWinThick = th; }
    void SetTargetLength(G4double l) {fTarLength = l; }
    void SetBeamlineOn(G4int b) {fBeamline = b; }

    G4OpticalSurface* GetSurface(void) { return fSurface; }

    void AddWorldMPV(const G4String& prop, G4MaterialPropertyVector* mpv);
    void AddWorldMPC(const G4String& prop, G4double v);
    G4MaterialPropertiesTable* GetWorldMaterialPropertiesTable()
    {
        return fWorldMPT;
    }

  void SetWorldMaterial(const G4String&);
  G4Material* GetWorldMaterial() const { return fWorldMaterial; }

private:
    G4NistManager* fNistManager;
    DetectorMessenger* fDetMessenger;
    
    G4double fNPSAngle;
    G4double fNPSDist;
    G4double fNPSShieldThick;
    G4double fHCALAngle;
    G4double fHCALDist;
    G4double fHCALShieldThick;
    G4double fSCWinThick;
    G4double fTarLength;
    G4double fBeamline;
    
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;
    
    static const G4int fNPSrow = 5;
    static const G4int fNPSNcol = 240;
    
    static const G4int fHCALNrow = 3;
    static const G4int fHCALNcol = 96;
    
    static const G4int fHodoNrow = 15;
    static const G4int fHodoNcol = 480;
    
    static const G4int fNSD = ( fNPSNrow*fNPSNcol + fHodoNrow*fHodoNcol + fHCALNrow*HCALNcol + 1);

    G4VPhysicalVolume* world_PV;
    G4VPhysicalVolume* fDetVol[fNSD];
    G4LogicalVolume* fWorld_LV;
    G4LogicalVolume* fLogicTarget;

    G4Material* fWorldMaterial;
    G4MaterialPropertiesTable* fWorldMPT;
    
    VirtualDetectorSD* fVirtualDetectorSD;
    RealDetectorSD* fRealDetectorSD;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
