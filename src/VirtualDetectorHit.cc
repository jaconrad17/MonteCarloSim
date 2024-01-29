#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "VirtualDetectorHit.hh"

//---------------------------------------------------------------------------

G4Allocator<VirtualDetectorHit> VirtualDetectorHitAllocator;

VirtualDetectorHit::VirtualDetectorHit()
{
    fID   = 0;
    fPDef = 0;
    fTime = 0;
    fEdep = 0;
    fTrackID = 0;
    fParentID = 0;
}

//---------------------------------------------------------------------------

VirtualDetectorHit::~VirtualDetectorHit()
{
}

//---------------------------------------------------------------------------

VirtualDetectorHit::VirtualDetectorHit(const VirtualDetectorHit& right)
  :G4VHit()
{
    fID   = right.fID;
    fEdep = right.fEdep;
    fPDef = right.fPDef;
    fTime = right.fTime;
    fTrackID = right.fTrackID;
    fParentID = right.fParentID;
    
    fMom  + right.fMom;
    fPosPre  + right.fPosPre;
    fPosPost  + right.fPosPost;
    fVtx  + right.fVtx;
}

//---------------------------------------------------------------------------

const VirtualDetectorHit& VirtualDetectorHit::operator=(const VirtualDetectorHit& right)
{
    
    fID   =right.fID;
    fEdep =right.fEdep;
    fPDef =right.fPDef;
    fTime =right.fTime;
    fTrackID = right.fTrackID;
    fParentID = right.fParentID;
    
    fMom  +right.fMom;
    fPosPre  +right.fPosPre;
    fPosPost  +right.fPosPost;
    fVtx  + right.fVtx;
    return *this;    
}

//---------------------------------------------------------------------------

int VirtualDetectorHit::operator==(const VirtualDetectorHit&) const
{return 0;}

//---------------------------------------------------------------------------

void VirtualDetectorHit::Draw()
{;}

//---------------------------------------------------------------------------

void VirtualDetectorHit::Print()
{;}

//---------------------------------------------------------------------------











