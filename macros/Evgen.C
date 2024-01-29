#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TF1.h"
#include "TH1.h"

#include "cteq/cteqpdf.h"
#include "Inelastic.h"

#include <iostream>
using namespace std;

// ----------------------------------------------------------------------------

void   InitOutput( TString );
void   GenerateReaction();

// ----------------------------------------------------------------------------

cteq_pdf_t *__dis_pdf;
void   initcteqpdf();
double dissigma( double, double, double );

// ----------------------------------------------------------------------------

const Int_t     fMaxparticles = 10;
TVector3*       fVertex[fMaxparticles];
TLorentzVector* fP4Lab[fMaxparticles];

TRandom3*     fRand;
TDatabasePDG* fDBpdg;
TFile*        fROOTFile;
TTree*        fROOTTree;
Int_t         fNparticles;  
Double_t      fWeight;
Int_t         fReactFlag;  
Double_t      fVx[fMaxparticles];
Double_t      fVy[fMaxparticles];
Double_t      fVz[fMaxparticles];
Double_t      fPx[fMaxparticles];
Double_t      fPy[fMaxparticles];
Double_t      fPz[fMaxparticles];
Double_t      fE[fMaxparticles];
Int_t         fPDG[fMaxparticles];

enum { kElastic, kQE, kInelastic, kPiPhoto, kDIS, kNReact }; 

Double_t fRasterSize  = 0.2;   // cm
Double_t fTarLength   = 10.;   // cm
Double_t fWindowThick = 0.125; // cm
Double_t fBeamE       = 6.6;   // GeV
Double_t fThetaMin    = 10.5 * TMath::DegToRad(); 
Double_t fThetaMax    = 20.5 * TMath::DegToRad(); 

// ----------------------------------------------------------------------------

void Evgen( Int_t run_no = 9999, Int_t ngen = 1000 )
{
    
    TString outrootfile = Form("out/gen_%d.root", run_no );
    
    initcteqpdf();
    
    fRand = new TRandom3(0);
    
    InitOutput( outrootfile );
    
    for( int i=0 ; i < fMaxparticles ; i++ ) {
        fVertex[i] = new TVector3( 0., 0., -30. );
        fP4Lab[i]  = new TLorentzVector( 0., 0., fBeamE, fBeamE );
        fPDG[i]    = 11;
    }
    
    fDBpdg = new TDatabasePDG();
    TString pdgtable = gSystem->Getenv("ROOTSYS");
    pdgtable.Append("/etc/pdg_table.txt");
    fDBpdg->ReadPDGTable(pdgtable);
    
    // ----------------------------------------------------------------------------
    
    for( Int_t i = 0; i < ngen; i++ ) {
        
        if( i%1000 == 0 ) {
            printf("Event %8d\r", i);
            fflush(stdout);
        }
        
        fReactFlag = fRand->Integer( kNReact );
        GenerateReaction();
        
        if( fNparticles < 1 ) continue;
        
        for( int i=0 ; i < fNparticles ; i++ ) {
            
            fVx[i]  = (Double_t)fVertex[i]->X();
            fVy[i]  = (Double_t)fVertex[i]->Y();
            fVz[i]  = (Double_t)fVertex[i]->Z();
            
            fPx[i]  = 1000*(Double_t)fP4Lab[i]->Px();
            fPy[i]  = 1000*(Double_t)fP4Lab[i]->Py();
            fPz[i]  = 1000*(Double_t)fP4Lab[i]->Pz();
            fE[i]   = 1000*(Double_t)fP4Lab[i]->E();
        }
        
        fROOTTree->Fill();
    }
    
    // ----------------------------------------------------------------------------
    
    fROOTTree->Write();
    fROOTFile->Close();
        
}

// ----------------------------------------------------------------------------

void InitOutput( TString outname )
{
    
    fROOTFile = new TFile(outname,"RECREATE");
    fROOTTree = new TTree("TGen", "Generator tree");
    fROOTTree->SetAutoSave();
    
    fROOTTree->Branch("weight",     &fWeight,     "weight/D");
    fROOTTree->Branch("flag",       &fReactFlag,  "flag/I");
    fROOTTree->Branch("Nparticles", &fNparticles, "Nparticles/I");
    
    fROOTTree->Branch("vx",  fVx,  "vx[Nparticles]/D");
    fROOTTree->Branch("vy",  fVy,  "vy[Nparticles]/D");
    fROOTTree->Branch("vz",  fVz,  "vz[Nparticles]/D");
    fROOTTree->Branch("px",  fPx,  "px[Nparticles]/D");
    fROOTTree->Branch("py",  fPy,  "py[Nparticles]/D");
    fROOTTree->Branch("pz",  fPz,  "pz[Nparticles]/D");
    fROOTTree->Branch("E",   fE,   "E[Nparticles]/D");
    fROOTTree->Branch("pdg", fPDG, "pdg[Nparticles]/I");
}

// ----------------------------------------------------------------------------

void GenerateReaction()
{
    
    Double_t xv = fRand->Uniform( -fRasterSize/2, fRasterSize/2 );
    Double_t yv = fRand->Uniform( -fRasterSize/2, fRasterSize/2 );
    Double_t zv = fRand->Uniform( -fTarLength/2, fTarLength/2 );
    
    Double_t L      = (0.0708*6.022e23)*fTarLength*(65.e-6/1.602e-19); // 65uA on LH2
    Double_t dOmega = TMath::TwoPi()*(TMath::Cos(fThetaMin)-TMath::Cos(fThetaMax));
    
    TLorentzVector beamP4( 0., 0., fBeamE, fBeamE );
    TLorentzVector targP4( 0., 0., 0., fDBpdg->GetParticle(2212)->Mass() );
    
    Double_t mt = fDBpdg->GetParticle(2212)->Mass();
    Double_t me = fDBpdg->GetParticle(11)->Mass();
    Double_t th, ph;
    
    TVector3 scatP3;
    Double_t scatE, scatP, Q2;
    
    TLorentzVector scatP4, qP4;
    
    Double_t hbarc = 197.327/1000;
    Double_t alpha = 1/137.036;
    
    // ----------------------------------------------------------------------------
    // ep -> ep Kelly fit
    // ----------------------------------------------------------------------------
    
    if( fReactFlag == kElastic ) {
        
        fNparticles = 2;
        fVertex[0]->SetXYZ( xv, yv, zv );
        fVertex[1]->SetXYZ( xv, yv, zv );
        fPDG[0] = 11;
        fPDG[1] = 2212;
        
        Double_t mp = fDBpdg->GetParticle(fPDG[1])->Mass();
        
        th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
        ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
        scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );
        
        scatE = (fBeamE*targP4.E())/(fBeamE*(1 -TMath::Cos(th)) + targP4.E()-targP4.Vect().Dot(scatP3));
        scatP = TMath::Sqrt( scatE*scatE - me*me );
        scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE );
        qP4   = beamP4 - scatP4;
        Q2    = -qP4.M2();
        
        TLorentzVector recoilP4 = targP4 + qP4;
        *fP4Lab[0] = scatP4;
        *fP4Lab[1] = recoilP4;
        
        Double_t tau   = Q2/(4*mp*mp);
        Double_t GE    = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
        Double_t GM    = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
        
        Double_t dSigMott  = 1e7*hbarc*hbarc*alpha*alpha/
        (4*beamP4.E()*beamP4.E()*TMath::Power(TMath::Sin(th/2),4))
        *(scatP4.E()/beamP4.E())*TMath::Power(TMath::Cos(th/2),2);
        Double_t dSigRosen = dSigMott *
        ( (GE*GE + tau*GM*GM)/(1+tau) + (2*tau*GM*GM)*TMath::Power(TMath::Tan(th/2),2) ); // nb/sr
        
        fWeight = L * dOmega * dSigRosen*1e-33;
        //    cout << "EP " << fWeight << "\t" << L << "\t"<< dOmega << "\t" << dSigRosen*1e-33 << endl;
    }
    
    
    // ----------------------------------------------------------------------------
    // ep -> e(X) DIS with F2 from CTEQ6
    // ----------------------------------------------------------------------------
    
    if ( fReactFlag == kDIS ) {
        
        fNparticles = 1;
        fVertex[0]->SetXYZ( xv, yv, zv );
        fPDG[0] = 11;
        
        th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
        ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
        scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );
        
        Double_t Emin = 0.2;
        Double_t Emax = fBeamE / (1.0 + fBeamE/mt*(1.0-TMath::Cos(th))) - 0.001;
        
        scatE = fRand->Uniform( Emin, Emax );
        scatP = TMath::Sqrt( scatE*scatE - me*me );
        scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE );
        
        *fP4Lab[0] = scatP4;
        
        Double_t  dSigDIS = dissigma( (double)beamP4.E(), (double)th, (double)scatP4.E() ); // nb/(GeV sr)
        
        fWeight = L * dOmega * (Emax - Emin) * dSigDIS * 1e-33;
        //      cout << "DIS " << fWeight << "\t" << L << "\t"<< dOmega << "\t" << Emax-Emin << "\t" << dSigDIS*1e-33  << endl;
    }
    
    // ----------------------------------------------------------------------------
    // eN -> eN (Quasi-elastic in the target windows, Kelly fit)
    // ----------------------------------------------------------------------------
    
    if( fReactFlag == kQE ) {
        
        L = (2.7*6.022e23/27)*2*fWindowThick*(60.e-6/1.602e-19); // 60uA on Al windows
        
        if( fRand->Uniform(0,1) < 1.0/2.0 )
            zv = fRand->Uniform( -fWindowThick - fTarLength/2, -fTarLength/2 );
        else
            zv = fRand->Uniform( fTarLength/2, fTarLength/2 + fWindowThick);
        
        fNparticles = 2;
        fVertex[0]->SetXYZ( xv, yv, zv );
        fVertex[1]->SetXYZ( xv, yv, zv );
        
        Int_t pdg_N;
        if( fRand->Uniform(0,1) < 1.0/2.0 )
            pdg_N = 2212; // p
        else
            pdg_N = 2112; // n
        fPDG[0] = 11;
        fPDG[1] = pdg_N;
        
        Double_t mt       = fDBpdg->GetParticle(fPDG[1])->Mass();
        Double_t mN       = mt;
        Double_t pF       = fRand->Uniform( 0.0, 0.3  );
        Double_t EF       = TMath::Sqrt(pF*pF + mt*mt );
        Double_t costhtar = fRand->Uniform( -1., 1. );
        Double_t thtar    = TMath::ACos( costhtar );
        Double_t phtar    = fRand->Uniform( -TMath::Pi(), TMath::Pi() );
        
        TLorentzVector targP4( pF * TMath::Sin( thtar ) * TMath::Cos( phtar ),
                              pF * TMath::Sin( thtar ) * TMath::Sin( phtar ),
                              pF * TMath::Cos( thtar ),
                              EF );
        
        th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
        ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
        scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );
        
        scatE = (fBeamE*targP4.E())/(fBeamE*(1 -TMath::Cos(th)) + targP4.E()-targP4.Vect().Dot(scatP3));
        scatP = TMath::Sqrt( scatE*scatE - me*me );
        scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE );
        qP4   = beamP4 - scatP4;
        Q2    = -qP4.M2();
        
        TLorentzVector recoilP4 = targP4 + qP4;
        *fP4Lab[0] = scatP4;
        *fP4Lab[1] = recoilP4;
        
        Double_t tau   = Q2/(4*mN*mN);
        Double_t GE, GM, GD;
        if( pdg_N == 2212 ) {
            GE = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
            GM = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
        }
        else {
            GD = pow(1.0 + Q2/(0.71), -2.0);
            GE = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
            GM = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );
        }
        
        Double_t dSigMott  = 1e10*hbarc*hbarc*alpha*alpha/
        (4*beamP4.E()*beamP4.E()*TMath::Power(TMath::Sin(th/2),4))
        *(scatP4.E()/beamP4.E())*TMath::Power(TMath::Cos(th/2),2);
        Double_t dSigRosen = dSigMott *
        ( (GE*GE + tau*GM*GM)/(1+tau) + (2*tau*GM*GM)*TMath::Power(TMath::Tan(th/2),2) );
        
        fWeight = L * dOmega * dSigRosen*1e-36;
        //    cout << "QE " << fWeight << "\t" << L << "\t"<< dOmega << "\t" << dSigRosen*1e-36 << endl;
    }
    
    // ----------------------------------------------------------------------------
    // ep -> eppi0 and ep -> enpi+ (Bosted and Christy, arxiv 0712.3731v4 )
    // ----------------------------------------------------------------------------
    
    if ( fReactFlag == kInelastic ) {
        
        fNparticles = 3;
        fVertex[0]->SetXYZ( xv, yv, zv );
        fVertex[1]->SetXYZ( xv, yv, zv );
        fVertex[2]->SetXYZ( xv, yv, zv );
        
        th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
        ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
        scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );
        
        Double_t W2min = pow(( mt + fDBpdg->GetParticle(111)->Mass() ), 2.0);
        Double_t Emin  = 0.2;
        Double_t Emax  = (fBeamE / (1.0 + fBeamE/mt*(1.0-TMath::Cos(th))))-0.14;
        
        scatE = fRand->Uniform( Emin, Emax );
        scatP = TMath::Sqrt( scatE*scatE - me*me );
        scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE );
        qP4   = beamP4 - scatP4;
        Q2    = -qP4.M2();
        
        TLorentzVector XP4 = targP4 + qP4;
        Double_t       W2  = XP4.M2();
        Double_t       W   = TMath::Sqrt( W2 );
        Double_t       x   = -qP4.M2()/(2.0*targP4.Dot( qP4 ) );
        
        Int_t pdg_Nf, pdg_pi;
        if( fRand->Uniform(0,1) < 2.0/3.0 ){
            pdg_Nf = 2212; // p
            pdg_pi = 111;  // pi0
        } else {
            pdg_Nf = 2112; // n
            pdg_pi = 211;  // pi+
        }
        fPDG[1] = pdg_Nf;
        fPDG[2] = pdg_pi;
        
        Double_t M_ni = fDBpdg->GetParticle(2212)->Mass();
        Double_t M_nf = fDBpdg->GetParticle(pdg_Nf)->Mass();
        Double_t Mpi  = fDBpdg->GetParticle(pdg_pi)->Mass();
        
        Double_t thpi = TMath::ACos( fRand->Uniform(-1,1) );
        Double_t phpi = fRand->Uniform(0, TMath::TwoPi());
        
        Float_t Epcm  = (W2 + Mpi*Mpi - M_nf*M_nf)/(2.0*W);
        Float_t ppcm  = TMath::Sqrt( (Epcm*Epcm - Mpi*Mpi) );
        Float_t pxpcm = ppcm * TMath::Sin( thpi ) * TMath::Cos( phpi );
        Float_t pypcm = ppcm * TMath::Sin( thpi ) * TMath::Sin( phpi );
        Float_t pzpcm = -ppcm * TMath::Cos( thpi );
        
        TLorentzVector pionP4( pxpcm,pypcm,pzpcm,Epcm );
        TLorentzVector nucleonP4( -pxpcm,-pypcm,-pzpcm, TMath::Sqrt( ppcm*ppcm + M_nf*M_nf ) );
        
        TVector3 boost = XP4.BoostVector();
        pionP4.Boost( boost );
        nucleonP4.Boost( boost );
        
        *fP4Lab[0] = scatP4;
        *fP4Lab[1] = nucleonP4;
        *fP4Lab[2] = pionP4;
        
        Double_t dSigInel = sigma_p( (double)beamP4.E(), (double)th, (double)scatP4.E() );  // nd/(GeV sr)
        
        fWeight = L * dOmega * (Emax - Emin) * dSigInel * 1e-33;
        //    cout << "Inelastic " << fWeight << "\t" << L << "\t"<< dOmega << "\t" << Emax-Emin << "\t" << dSigInel*1e-33  << endl;
    }
    
    // ----------------------------------------------------------------------------
    // gammap -> ppi0 (fit to E99114 data)
    // ----------------------------------------------------------------------------
    
    if ( fReactFlag == kPiPhoto ) {
        
        fNparticles = 2;
        fVertex[0]->SetXYZ( xv, yv, zv );
        fVertex[1]->SetXYZ( xv, yv, zv );
        fPDG[0] = 2212;
        fPDG[1] = 111;
        
        th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
        ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
        
        Double_t Emin  = 2.0;
        Double_t Emax  = fBeamE;
        TF1* fbrem     = new TF1("fbrem","1/x", Emin, Emax);
        TH1* hbrem     = (TH1*)fbrem->GetHistogram();
        Double_t beamE = hbrem->GetRandom(); //fRand->Uniform( Emin, Emax );
        beamP4.SetPxPyPzE( 0., 0., beamE, beamE );
        
        TLorentzVector XP4 = targP4 + beamP4;
        Double_t       W2  = XP4.M2();
        Double_t       W   = TMath::Sqrt( W2 );
        
        Double_t M_ni = fDBpdg->GetParticle(2212)->Mass();
        Double_t M_nf = fDBpdg->GetParticle(fPDG[0])->Mass();
        Double_t Mpi  = fDBpdg->GetParticle(fPDG[1])->Mass();
        
        Double_t thp  = TMath::ACos( fRand->Uniform(0.47,0.67) );
        Double_t php  = fRand->Uniform(0, TMath::TwoPi());
        Float_t Epcm  = (W2 + M_nf*M_nf - Mpi*Mpi)/(2.0*W);
        Float_t ppcm  = TMath::Sqrt( (Epcm*Epcm - M_nf*M_nf) );
        Float_t pxpcm = ppcm * TMath::Sin( thp ) * TMath::Cos( php );
        Float_t pypcm = ppcm * TMath::Sin( thp ) * TMath::Sin( php );
        Float_t pzpcm = -ppcm * TMath::Cos( thp );
        
        TLorentzVector nucleonP4( pxpcm,pypcm,pzpcm,Epcm );
        TLorentzVector pionP4( -pxpcm,-pypcm,-pzpcm, TMath::Sqrt( ppcm*ppcm + Mpi*Mpi ) );
        
        TVector3 boost = XP4.BoostVector();
        pionP4.Boost( boost );
        nucleonP4.Boost( boost );
        
        *fP4Lab[0] = nucleonP4;
        *fP4Lab[1] = pionP4;
        
        Double_t costhpcm  = TMath::Cos( thp );
        Double_t conv      = TMath::Pi()/( nucleonP4.E() * nucleonP4.E() );
        Double_t Mandel_s  = 2*mt*beamE + mt*mt;
        
        Double_t dSigPiPhoto = conv * 0.72 * TMath::Power( 10.9/Mandel_s, 7.5 )
        * TMath::Power((1-costhpcm), -1.0); // nb/sr
        
        fWeight = L * 0.03 * dOmega * dSigPiPhoto * 1e-33; // need brem photon frac correction
        //    cout << "Pion photo " << fWeight << "\t" << L << "\t"<< dOmega << "\t" << dSigPiPhoto*1e-36 << endl;
        
    }
}

// ----------------------------------------------------------------------------
// Use CTEQ6 parameterization for DIS PDFs
// ----------------------------------------------------------------------------

void initcteqpdf(){
    
    __dis_pdf = cteq_pdf_alloc_id(400);
    
    assert(__dis_pdf);    
}

double dissigma( double ebeam, double th, double eprime ){
    
    // Return in nb/(GeV*sr)
    
    double Q2 = 2.0*eprime*ebeam*(1.0-cos(th));
    double nu = ebeam-eprime;
    double Mp = 0.93827;
    
    double x = Q2/(2.0*Mp*nu);
    double y = nu/ebeam;
    
    double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
    double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
    double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
    double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );
    
    double quv = qu-qubar;
    double qdv = qd-qdbar;
    
    double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );
    
    double F2 = 0.0;
    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;
    
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv );
    F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
    
    double F1 = F2/(2.0*x);
    
    double hbarc = 197.327/1000;
    double alpha = 1/137.036;
    
    // From PDG
    double ds_dxdy = 4.0*TMath::Pi()*pow(alpha,2)*((1.0-y-pow(x*y*Mp,2)/Q2)*F2+y*y*x*F1)/(x*y*Q2);
    
    // In GeV^-2
    double ds_dOmega_dE = ds_dxdy*eprime/(2.0*TMath::Pi()*Mp*nu);
    
    return ds_dOmega_dE*pow(hbarc,2)*1e7; // GeV2 -> nb
}

// ----------------------------------------------------------------------------
