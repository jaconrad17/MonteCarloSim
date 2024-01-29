#ifndef OutputManager_h
#define OutputManager_h 1

#include "globals.hh"
#include "PrimaryGeneratorAction.hh"
#include "OutputMessenger.hh"

#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "Rtypes.h"
#include "TVector3.h"
#include "TString.h"

class TTree;
class TFile;

//---------------------------------------------------------------------------

class OutputManager {
    
public:
    
    OutputManager();
    ~OutputManager();
    
    void InitOutput();
    
    void ZeroArray();
    void FillPrimaryArray( Int_t );
    void FillVirtualArray( Int_t );
    void FillRealArray( Int_t );
    void FillTree( G4double, G4int, G4int );
    
    void SetOutFileName( TString fname ) { fOutFileName  = fname; }
    
    void SetPrimaryEnergy   (G4double       ene  )       { fPrimary_energy = ene;  }
    void SetPrimaryPDef     (G4ParticleDefinition* pdef) { fPrimary_PDef   = pdef; }
    void SetPrimaryDirection(G4ThreeVector  dir  )       { fPrimary_dir    = dir;  }
    void SetPrimaryVertex   (G4ThreeVector  vtx  )       { fPrimary_vtx    = vtx;  }
    
    void SetVirtualPDef      ( G4ParticleDefinition* sp ) { fVirtual_pdef = sp;    }
    void SetVirtualPosPre    ( G4ThreeVector  spos )      { fVirtual_pospre  = spos;  }
    void SetVirtualPosPost   ( G4ThreeVector  spos )      { fVirtual_pospost  = spos;  }
    void SetVirtualVertex    ( G4ThreeVector  vtx )       { fVirtual_vtx  = vtx;  }
    void SetVirtualP3        ( G4ThreeVector  smom )      { fVirtual_p3   = smom;  }
    void SetVirtualTime      ( G4double       stime )     { fVirtual_time = stime; }
    void SetVirtualID        ( G4int          sid )       { fVirtual_detid   = sid;   }
    void SetVirtualEnergy    ( G4double       ene)        { fVirtual_energy = ene; }
    void SetVirtualPID       ( G4int          pid )       { fVirtual_Pid   = pid;   }
    void SetVirtualTID       ( G4int          tid )       { fVirtual_Tid   = tid;   }
    
    void SetRealPosPre    ( G4ThreeVector  epos )      { fReal_pospre  = epos;  }
    void SetRealPosPost   ( G4ThreeVector  epos )      { fReal_pospost  = epos;  }
    void SetRealID        ( G4int          eid )       { fReal_detid   = eid;   }
    void SetRealEdep      ( G4double       eedep)      { fReal_edep = eedep; }
    void SetRealTime      ( G4double       time)       { fReal_time = time; }
    
private:
    
    OutputMessenger*      fOutMessenger;
    TString               fOutFileName;
    TFile*                fROOTfile;
    TTree*                fROOTtree;
    
    const G4int           fNearm    = 1200;
    const G4int           fNearmCol = 240;
    const G4int           fNhodo    = 7200;
    const G4int           fNhodoCol = 480;
    const G4int           fNharm    = 288;
    const G4int           fNharmCol  = 96;
    
    // Event
    Double_t               fEvent_weight;
    Int_t                 fEvent_flag;
    Int_t                 fEvent_num;
    
    // Primary
    static const Int_t    fMaxprim = 50;
    
    G4ParticleDefinition* fPrimary_PDef;
    G4ThreeVector         fPrimary_dir;
    G4ThreeVector         fPrimary_vtx;
    G4double              fPrimary_energy;
    
    Int_t                 fPrimary_Nhits;
    Int_t                 fPrimary_pdg[fMaxprim];
    Double_t               fPrimary_E[fMaxprim];
    Double_t               fPrimary_xpos[fMaxprim];
    Double_t               fPrimary_ypos[fMaxprim];
    Double_t               fPrimary_zpos[fMaxprim];
    Double_t               fPrimary_px[fMaxprim];
    Double_t               fPrimary_py[fMaxprim];
    Double_t               fPrimary_pz[fMaxprim];
    
    // Virtual detectors
    static const Int_t    fMaxhits = 50000;
    
    G4ParticleDefinition* fVirtual_pdef;
    G4ThreeVector         fVirtual_p3;
    G4ThreeVector         fVirtual_pospre;
    G4ThreeVector         fVirtual_pospost;
    G4ThreeVector         fVirtual_vtx;
    G4double              fVirtual_time;
    G4int                 fVirtual_detid;
    G4int                 fVirtual_Pid;
    G4int                 fVirtual_Tid;
    G4double              fVirtual_energy;
    
    Int_t                 fVirtual_Nhits;
    Int_t                 fVirtual_pdg[fMaxhits];
    Int_t                 fVirtual_det[fMaxhits];
    Int_t                 fVirtual_row[fMaxhits];
    Int_t                 fVirtual_col[fMaxhits];
    Int_t                 fVirtual_tid[fMaxhits];
    Int_t                 fVirtual_pid[fMaxhits];
    Double_t               fVirtual_E[fMaxhits];
    Double_t               fVirtual_t[fMaxhits];
    Double_t               fVirtual_xpos[fMaxhits];
    Double_t               fVirtual_ypos[fMaxhits];
    Double_t               fVirtual_zpos[fMaxhits];
    Double_t               fVirtual_px[fMaxhits];
    Double_t               fVirtual_py[fMaxhits];
    Double_t               fVirtual_pz[fMaxhits];
    Double_t               fVirtual_vx[fMaxhits];
    Double_t               fVirtual_vy[fMaxhits];
    Double_t               fVirtual_vz[fMaxhits];
    
    // Real detectors
    G4ThreeVector         fReal_pospre;
    G4ThreeVector         fReal_pospost;
    G4int                 fReal_detid;
    G4double              fReal_edep;
    G4double              fReal_time;
    
    Int_t                 fReal_Nhits;
    Int_t                 fReal_det[fMaxhits];
    Int_t                 fReal_row[fMaxhits];
    Int_t                 fReal_col[fMaxhits];
    Double_t               fReal_Edep[fMaxhits];
    Double_t               fReal_t[fMaxhits];
    Double_t               fReal_xpos[fMaxhits];
    Double_t               fReal_ypos[fMaxhits];
    Double_t               fReal_zpos[fMaxhits];
    
};

#endif

//---------------------------------------------------------------------------
