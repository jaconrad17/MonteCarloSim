#ifndef VirtualDetectorHit_h
#define VirtualDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

//---------------------------------------------------------------------------

class VirtualDetectorHit : public G4VHit
{
public:
    
    VirtualDetectorHit();
    ~VirtualDetectorHit();
    VirtualDetectorHit(const VirtualDetectorHit&);
    const VirtualDetectorHit& operator=(const
                                        VirtualDetectorHit&);
    int operator==(const VirtualDetectorHit&) const;
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
    virtual void Draw();
    void Print();
    
protected:
    
    G4int                  fID;
    G4ParticleDefinition*  fPDef;
    G4double               fTime;
    G4double               fEdep;
    G4ThreeVector          fMom;
    G4ThreeVector          fPosPre;
    G4ThreeVector          fPosPost;
    G4ThreeVector          fVtx;
    G4int                  fTrackID;
    G4int                  fParentID;
    
public:
    
    inline void SetID       (G4int i)                  { fID  = i;   };
    inline void SetPDef     (G4ParticleDefinition* pd) { fPDef = pd;  };
    inline void SetTime     (G4double t)               { fTime = t;   };
    inline void SetMomentum (G4ThreeVector p)          { fMom = p;    };
    inline void SetPrePosition (G4ThreeVector pos)     { fPosPre=pos;    };
    inline void SetPostPosition (G4ThreeVector pos)    { fPosPost=pos;    };
    inline void SetVertex   (G4ThreeVector vtx)        { fVtx=vtx;    };
    inline void SetEnergy   (G4double de)              { fEdep = de; };
    inline void SetTrackID    (G4int tt)               { fTrackID = tt; };
    inline void SetParentID   (G4int pp)               { fParentID = pp; };
    
    inline G4int                 GetID()               { return fID;  };
    inline G4int                 GetTrackID()          { return fTrackID;  };
    inline G4int                 GetParentID()         { return fParentID;  };
    inline G4ParticleDefinition* GetPDef()             { return fPDef; };
    inline G4double              GetTime()             { return fTime; };
    inline G4ThreeVector         GetMom()              { return fMom; };
    inline G4ThreeVector         GetPosPre()           { return fPosPre; };
    inline G4ThreeVector         GetPosPost()          { return fPosPost; };
    inline G4ThreeVector         GetVertex()           { return fVtx; };
    inline G4double              GetEdep()             { return fEdep; };
};

//---------------------------------------------------------------------------

typedef G4THitsCollection<VirtualDetectorHit> VirtualDetectorHitsCollection;

extern G4Allocator<VirtualDetectorHit> VirtualDetectorHitAllocator;


inline void* VirtualDetectorHit::operator new(size_t)
{
    void* aHit;
    aHit = (void*) VirtualDetectorHitAllocator.MallocSingle();
    return aHit;
}


inline void VirtualDetectorHit::operator delete(void* aHit)
{
    VirtualDetectorHitAllocator.FreeSingle((VirtualDetectorHit*) aHit);
}

#endif

//---------------------------------------------------------------------------










