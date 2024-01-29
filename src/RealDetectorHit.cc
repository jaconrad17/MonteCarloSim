#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "RealDetectorHit.hh"

//---------------------------------------------------------------------------

G4Allocator<RealDetectorHit> RealDetectorHitAllocator;

RealDetectorHit::RealDetectorHit()
{
    fID   = 0;
    fEdep = 0;
    fTime = 0;
}

//---------------------------------------------------------------------------

RealDetectorHit::~RealDetectorHit()
{
}

//---------------------------------------------------------------------------

RealDetectorHit::RealDetectorHit(const RealDetectorHit& right)
  :G4VHit()
{
    fID   = right.fID;
    fEdep = right.fEdep;
    fTime = right.fTime;
    
    fPosPre  + right.fPosPre;
    fPosPost  + right.fPosPost;
}

//---------------------------------------------------------------------------

const RealDetectorHit& RealDetectorHit::operator=(const RealDetectorHit& right)
{
    
    fID   =right.fID;
    fEdep =right.fEdep;
    fTime = right.fTime;
    
    fPosPre  +right.fPosPre;
    fPosPost  +right.fPosPost;
    return *this;
}

//---------------------------------------------------------------------------

int RealDetectorHit::operator==(const RealDetectorHit&) const
{return 0;}

//---------------------------------------------------------------------------

void RealDetectorHit::Draw()
{;}

//---------------------------------------------------------------------------

void RealDetectorHit::Print()
{;}

//---------------------------------------------------------------------------











