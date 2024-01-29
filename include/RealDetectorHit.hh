#ifndef RealDetectorHit_h
#define RealDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"

//---------------------------------------------------------------------------

class RealDetectorHit : public G4VHit
{
public:
    
    RealDetectorHit();
    ~RealDetectorHit();
    RealDetectorHit(const RealDetectorHit&);
    const RealDetectorHit& operator=(const
                                     RealDetectorHit&);
    int operator==(const RealDetectorHit&) const;
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
    virtual void Draw();
    void Print();
    
protected:
    
    G4int                  fID;
    G4double               fEdep;
    G4double               fTime;
    G4ThreeVector          fPosPre;
    G4ThreeVector          fPosPost;
    
public:
    
    inline void SetID       (G4int i)                  { fID  = i;   };
    inline void SetPrePosition (G4ThreeVector pos)     { fPosPre=pos;    };
    inline void SetPostPosition (G4ThreeVector pos)    { fPosPost=pos;    };
    inline void AddEnergy   (G4double de)              { fEdep += de; };
    inline void SetTime(G4double t) {fTime = t;};
    
    inline G4int                 GetID()               { return fID;  };
    inline G4ThreeVector         GetPosPre()           { return fPosPre; };
    inline G4ThreeVector         GetPosPost()          { return fPosPost; };
    inline G4double              GetEdep()             { return fEdep; };
    inline G4double    GetTime()   { return fTime; };
        
};

//---------------------------------------------------------------------------

typedef G4THitsCollection<RealDetectorHit> RealDetectorHitsCollection;

extern G4Allocator<RealDetectorHit> RealDetectorHitAllocator;


inline void* RealDetectorHit::operator new(size_t)
{
    void* aHit;
    aHit = (void*) RealDetectorHitAllocator.MallocSingle();
    return aHit;
}


inline void RealDetectorHit::operator delete(void* aHit)
{
    RealDetectorHitAllocator.FreeSingle((RealDetectorHit*) aHit);
}

#endif

//---------------------------------------------------------------------------










