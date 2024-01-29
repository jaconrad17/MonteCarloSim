#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

using namespace std;

const Bool_t ApplyThresh     = true;
const Bool_t ApplyWindow     = true;

const Double_t Earm_threshold = 0.5;
const Double_t Earm_window    = 5000;
const Double_t Harm_threshold = 0.5;
const Double_t Harm_window    = 5000;
const Double_t Hodo_threshold = 0.5;
const Double_t Hodo_window    = 5000;

const Int_t   NEarm          = 1201;
const Int_t   NHodo          = 7201;
const Int_t   NHarm          = 289;

void SinglesBG_OpNovice2( )
{
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TRandom3* fRand = new TRandom3(0);
    
    TChain* TOut = new TChain("TOut");
    
    TOut->Add("out/batch_990*.root"); // beam
    TOut->Add("out/batch_13000.root"); // elastic and DIS
    
    const Int_t     Maxprim = 50;
    const Int_t     Maxhits = 10000;
    
    Double_t         Event_weight;
    Int_t           Event_flag;
    Int_t           Event_num;
    Int_t           Primary_Nhits;
    Int_t           Primary_pdg[Maxprim];
    Double_t         Primary_E[Maxprim];
    Double_t         Primary_x[Maxprim];
    Double_t         Primary_y[Maxprim];
    Double_t         Primary_z[Maxprim];
    Double_t         Primary_px[Maxprim];
    Double_t         Primary_py[Maxprim];
    Double_t         Primary_pz[Maxprim];
    Int_t           Virtual_Nhits;
    Int_t           Virtual_pdg[Maxhits];
    Int_t           Virtual_det[Maxhits];
    Int_t           Virtual_row[Maxhits];
    Int_t           Virtual_col[Maxhits];
    Int_t           Virtual_tid[Maxhits];
    Int_t           Virtual_pid[Maxhits];
    Double_t         Virtual_E[Maxhits];
    Double_t         Virtual_t[Maxhits];
    Double_t         Virtual_x[Maxhits];
    Double_t         Virtual_y[Maxhits];
    Double_t         Virtual_z[Maxhits];
    Double_t         Virtual_px[Maxhits];
    Double_t         Virtual_py[Maxhits];
    Double_t         Virtual_pz[Maxhits];
    Double_t         Virtual_vx[Maxhits];
    Double_t         Virtual_vy[Maxhits];
    Double_t         Virtual_vz[Maxhits];
    Int_t           Real_Nhits;
    Int_t           Real_det[Maxhits];
    Int_t           Real_row[Maxhits];
    Int_t           Real_col[Maxhits];
    Double_t         Real_edep[Maxhits];
    Double_t         Real_t[Maxhits];
    Double_t         Real_x[Maxhits];
    Double_t         Real_y[Maxhits];
    Double_t         Real_z[Maxhits];
    
    TOut->SetBranchAddress("Event_weight",&Event_weight);
    TOut->SetBranchAddress("Event_flag",&Event_flag);
    TOut->SetBranchAddress("Event_num",&Event_num);
    TOut->SetBranchAddress("Primary_Nhits",&Primary_Nhits);
    TOut->SetBranchAddress("Primary_pdg",Primary_pdg);
    TOut->SetBranchAddress("Primary_E",Primary_E);
    TOut->SetBranchAddress("Primary_x",Primary_x);
    TOut->SetBranchAddress("Primary_y",Primary_y);
    TOut->SetBranchAddress("Primary_z",Primary_z);
    TOut->SetBranchAddress("Primary_px",Primary_px);
    TOut->SetBranchAddress("Primary_py",Primary_py);
    TOut->SetBranchAddress("Primary_pz",Primary_pz);
    TOut->SetBranchAddress("Virtual_Nhits",&Virtual_Nhits);
    TOut->SetBranchAddress("Virtual_pdg",Virtual_pdg);
    TOut->SetBranchAddress("Virtual_det",Virtual_det);
    TOut->SetBranchAddress("Virtual_row",Virtual_row);
    TOut->SetBranchAddress("Virtual_col",Virtual_col);
    TOut->SetBranchAddress("Virtual_tid",Virtual_tid);
    TOut->SetBranchAddress("Virtual_pid",Virtual_pid);
    TOut->SetBranchAddress("Virtual_E",Virtual_E);
    TOut->SetBranchAddress("Virtual_t",Virtual_t);
    TOut->SetBranchAddress("Virtual_x",Virtual_x);
    TOut->SetBranchAddress("Virtual_y",Virtual_y);
    TOut->SetBranchAddress("Virtual_z",Virtual_z);
    TOut->SetBranchAddress("Virtual_px",Virtual_px);
    TOut->SetBranchAddress("Virtual_py",Virtual_py);
    TOut->SetBranchAddress("Virtual_pz",Virtual_pz);
    TOut->SetBranchAddress("Virtual_vx",Virtual_vx);
    TOut->SetBranchAddress("Virtual_vy",Virtual_vy);
    TOut->SetBranchAddress("Virtual_vz",Virtual_vz);
    TOut->SetBranchAddress("Real_Nhits",&Real_Nhits);
    TOut->SetBranchAddress("Real_det",Real_det);
    TOut->SetBranchAddress("Real_row",Real_row);
    TOut->SetBranchAddress("Real_col",Real_col);
    TOut->SetBranchAddress("Real_edep",Real_edep);
    TOut->SetBranchAddress("Real_t",Real_t);
    TOut->SetBranchAddress("Real_x",Real_x);
    TOut->SetBranchAddress("Real_y",Real_y);
    TOut->SetBranchAddress("Real_z",Real_z);
    
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    
    gStyle->SetPadTopMargin(.05);
    gStyle->SetPadLeftMargin(.18);
    gStyle->SetPadRightMargin(.18);
    gStyle->SetPadBottomMargin(.15);
    
    gStyle->SetTitleOffset(1.1, "X");
    gStyle->SetTitleOffset(1.4, "Y");
    gStyle->SetTitleFont(42,"X");
    gStyle->SetTitleFont(42,"Y");
    gStyle->SetTitleSize(0.055,"X");
    gStyle->SetTitleSize(0.055,"Y");
    
    gStyle->SetLabelOffset(0.01, "X");
    gStyle->SetLabelOffset(0.01, "Y");
    gStyle->SetLabelFont(42,"X");
    gStyle->SetLabelFont(42,"Y");
    gStyle->SetLabelSize(0.045,"X");
    gStyle->SetLabelSize(0.045,"Y");
    
    gStyle->SetNdivisions(105,"X");
    gStyle->SetNdivisions(105,"Y");
    
    gStyle->SetStripDecimals(kFALSE);
    
    Long64_t nentries = TOut->GetEntries();
    Int_t ntrees      = TOut->GetNtrees();
    
    cout << "Processing TOut chain with " << ntrees << " trees and " << nentries << " total events" << endl;
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TH1F* hRealEarm_N   = new TH1F("hRealEarm_N",  "", NEarm, 0.,NEarm);
    TH1F* hRealEarm_Nb  = new TH1F("hRealEarm_Nb",  "", NEarm, 0.,NEarm);
    TH1F* hRealEarm_Ne  = new TH1F("hRealEarm_Ne",  "", NEarm, 0.,NEarm);
    TH1F* hRealEarm_Nd  = new TH1F("hRealEarm_Nd",  "", NEarm, 0.,NEarm);
    TH1F* hRealEarm_E   = new TH1F("hRealEarm_E",  "", 100,0.,5000.);
    TH1F* hRealEarm_Eb  = new TH1F("hRealEarm_Eb",  "", 100,0.,5000.);
    TH1F* hRealEarm_Ee  = new TH1F("hRealEarm_Ee",  "", 100,0.,5000.);
    TH1F* hRealEarm_Ed  = new TH1F("hRealEarm_Ed",  "", 100,0.,5000.);
    TH1F* hRealEarm_z   = new TH1F("hRealEarm_z",  "", 100,0.,800.);
    TH1F* hRealEarm_t   = new TH1F("hRealEarm_t",  "", 100, 0.,50.);
    TH2F* hRealEarm_xy  = new TH2F("hRealEarm_xy", "", 100,-500.,500., 100,-500.,500. );
    TH1F* hRealEarm_Edet[NEarm];
    for(Int_t i=0; i<NEarm; i++)
        hRealEarm_Edet[i] = new TH1F( Form("hRealEarm_Edet%d", i), "", 100,0.,1000.);
    
    TH1F* hRealHarm_N   = new TH1F("hRealHarm_N",  "", NHarm,0.,NHarm);
    TH1F* hRealHarm_Nb  = new TH1F("hRealHarm_Nb",  "", NHarm,0.,NHarm);
    TH1F* hRealHarm_Ne  = new TH1F("hRealHarm_Ne",  "", NHarm,0.,NHarm);
    TH1F* hRealHarm_Nd  = new TH1F("hRealHarm_Nd",  "", NHarm,0.,NHarm);
    TH1F* hRealHarm_E   = new TH1F("hRealHarm_E",  "", 100,0.,300.);
    TH1F* hRealHarm_Eb  = new TH1F("hRealHarm_Eb",  "", 100,0.,300.);
    TH1F* hRealHarm_Ee  = new TH1F("hRealHarm_Ee",  "", 100,0.,300.);
    TH1F* hRealHarm_Ed  = new TH1F("hRealHarm_Ed",  "", 100,0.,300.);
    TH1F* hRealHarm_z   = new TH1F("hRealHarm_z",  "", 100,0.,500.);
    TH1F* hRealHarm_t   = new TH1F("hRealHarm_t",  "", 100, 0.,50.);
    TH2F* hRealHarm_xy  = new TH2F("hRealHarm_xy", "", 100,-500.,500., 100,-500.,500. );
    TH1F* hRealHarm_Edet[NHarm];
    for(Int_t i=0; i<NHarm; i++)
        hRealHarm_Edet[i] = new TH1F( Form("hRealHarm_Edet%d", i), "", 100,0.,1000.);
    
    TH1F* hRealHodo_N   = new TH1F("hRealHodo_N",  "", 7200,0.,NHodo);
    TH1F* hRealHodo_Ne  = new TH1F("hRealHodo_Ne",  "", 7200,0.,NHodo);
    TH1F* hRealHodo_Nb  = new TH1F("hRealHodo_Nb",  "", 7200,0.,NHodo);
    TH1F* hRealHodo_Nd  = new TH1F("hRealHodo_Nd",  "", 7200,0.,NHodo);
    TH1F* hRealHodo_E   = new TH1F("hRealHodo_E",  "", 100,0.,50.0);
    TH1F* hRealHodo_Ee  = new TH1F("hRealHodo_Ee",  "", 100,0.,50.0);
    TH1F* hRealHodo_Eb  = new TH1F("hRealHodo_Eb",  "", 100,0.,50.0);
    TH1F* hRealHodo_Ed  = new TH1F("hRealHodo_Ed",  "", 100,0.,50.0);
    TH1F* hRealHodo_z   = new TH1F("hRealHodo_z",  "", 100,0.,500.);
    TH1F* hRealHodo_t   = new TH1F("hRealHodo_t",  "", 100, 0.,50.);
    TH2F* hRealHodo_xy  = new TH2F("hRealHodo_xy", "", 100,-500.,500., 100,-500.,500. );
    TH1F* hRealHodo_Edet[NHodo];
    for(Int_t i=0; i<NHodo; i++)
        hRealHodo_Edet[i] = new TH1F( Form("hRealHodo_Edet%d", i), "", 100,0.,0.3);
    
    TH1F* hEarm_E   = new TH1F("hEarm_E",  "", 100,3000.,6000.);
    TH1F* hEarm_Ee  = new TH1F("hEarm_Ee",  "", 100,3000.,6000.);
    TH1F* hEarm_Eq  = new TH1F("hEarm_Eq",  "", 100,3000.,6000.);
    TH1F* hEarm_Ei  = new TH1F("hEarm_Ei",  "", 100,3000.,6000.);
    TH1F* hEarm_Ep  = new TH1F("hEarm_Ep",  "", 100,3000.,6000.);
    
    TH1F* hEarmC_E   = new TH1F("hEarmC_E",  "", 100,3000.,6000.);
    TH1F* hEarmC_Ee  = new TH1F("hEarmC_Ee",  "", 100,3000.,6000.);
    TH1F* hEarmC_Eq  = new TH1F("hEarmC_Eq",  "", 100,3000.,6000.);
    TH1F* hEarmC_Ei  = new TH1F("hEarmC_Ei",  "", 100,3000.,6000.);
    TH1F* hEarmC_Ed  = new TH1F("hEarmC_Ed",  "", 100,3000.,6000.);
    TH1F* hEarmC_Ep  = new TH1F("hEarmC_Ep",  "", 100,3000.,6000.);
    
    TH1F* hEarmC4_E   = new TH1F("hEarmC4_E",  "", 100,4500.,6000.);
    TH1F* hEarmC4_Ee  = new TH1F("hEarmC4_Ee",  "", 100,4500.,6000.);
    TH1F* hEarmC4_Eq  = new TH1F("hEarmC4_Eq",  "", 100,4500.,6000.);
    TH1F* hEarmC4_Ei  = new TH1F("hEarmC4_Ei",  "", 100,4500.,6000.);
    TH1F* hEarmC4_Ed  = new TH1F("hEarmC4_Ed",  "", 100,4500.,6000.);
    TH1F* hEarmC4_Ep  = new TH1F("hEarmC4_Ep",  "", 100,4500.,6000.);
    
    TH1F* hEarmC5_E   = new TH1F("hEarmC5_E",  "", 100,5000.,6000.);
    TH1F* hEarmC5_Ee  = new TH1F("hEarmC5_Ee",  "", 100,5000.,6000.);
    TH1F* hEarmC5_Eq  = new TH1F("hEarmC5_Eq",  "", 100,5000.,6000.);
    TH1F* hEarmC5_Ei  = new TH1F("hEarmC5_Ei",  "", 100,5000.,6000.);
    TH1F* hEarmC5_Ed  = new TH1F("hEarmC5_Ed",  "", 100,5000.,6000.);
    TH1F* hEarmC5_Ep  = new TH1F("hEarmC5_Ep",  "", 100,5000.,6000.);
    
    TH1F* hEarmC_PhDiff   = new TH1F("hEarmC_PhDiff",  "",  100,-10.,10.);
    TH1F* hEarmC_PhDiffe  = new TH1F("hEarmC_PhDiffe",  "", 100,-10.,10.);
    TH1F* hEarmC_PhDiffq  = new TH1F("hEarmC_PhDiffq",  "", 100,-10.,10.);
    TH1F* hEarmC_PhDiffi  = new TH1F("hEarmC_PhDiffi",  "", 100,-10.,10.);
    TH1F* hEarmC_PhDiffd  = new TH1F("hEarmC_PhDiffd",  "", 100,-10.,10.);
    TH1F* hEarmC_PhDiffp  = new TH1F("hEarmC_PhDiffp",  "", 100,-10.,10.);
    
    TH1F* hEarmC_tDiff   = new TH1F("hEarmC_tDiff",  "",  100,-10.,10.);
    TH1F* hEarmC_tDiffe  = new TH1F("hEarmC_tDiffe",  "", 100,-10.,10.);
    TH1F* hEarmC_tDiffq  = new TH1F("hEarmC_tDiffq",  "", 100,-10.,10.);
    TH1F* hEarmC_tDiffi  = new TH1F("hEarmC_tDiffi",  "", 100,-10.,10.);
    TH1F* hEarmC_tDiffd  = new TH1F("hEarmC_tDiffd",  "", 100,-10.,10.);
    TH1F* hEarmC_tDiffp  = new TH1F("hEarmC_tDiffp",  "", 100,-10.,10.);
    
    TH1F* hEarmC4_Q2   = new TH1F("hEarmC4_Q2",  "",  100, 1.8, 3.2);
    TH1F* hEarmC5_Q2   = new TH1F("hEarmC5_Q2",  "",  100, 1.8, 3.2);
    
    TH1F* hEarmC4_Th   = new TH1F("hEarmC4_Th",  "",  100, 12.5, 18.5);
    TH1F* hEarmC5_Th   = new TH1F("hEarmC5_Th",  "",  100, 12.5, 18.5);
    
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    Int_t nhith  = 0;
    Int_t ngoodh   = 0;
    Int_t ngoodhod = 0;
    
    TLorentzVector p4h(0,0,0,0.938);
    TLorentzVector p4e(0,0,0,0.0);
    TLorentzVector Tp4(0,0,0,0.938); //target 4vec
    TLorentzVector kp4(0,0,6.6,6.6); //beam 4vec
    TLorentzVector Qp4, kpp4, Rp4; //q, recoil electron, recoil nucleon
    
    double timee, timeh;
    
    for(Long64_t ev=0; ev<nentries;ev++) {
        
        TOut->GetEntry(ev);
        
        if( ev%10000 == 0 ) {
            printf("Event %8lld\r", ev);
            fflush(stdout);
        }
        
        if( Event_flag == 999 )
            Event_weight = Event_weight * (1./(50.*Event_num));
        else
            Event_weight = Event_weight * (5./(50.*Event_num));
        
        if (Event_weight == 0 )
            Event_weight = 1;
        
        bool goodVe = false;
        bool goodVh = false;
        
        for( int i =0; i < Virtual_Nhits; i++ ) {
            
            double mom = TMath::Sqrt( Virtual_px[i]*Virtual_px[i] + Virtual_py[i]*Virtual_py[i] + Virtual_pz[i]*Virtual_pz[i] );
            
            if( Virtual_det[i] == 1 && Virtual_E[i] > 100 && Virtual_pdg[i] == 2212 ) {
                
                goodVh = true;
                
                p4h.SetPxPyPzE( 0.001*Virtual_px[i], 0.001*Virtual_py[i], 0.001*Virtual_pz[i], (0.001*Virtual_E[i] +0.938));
                
                timeh = Virtual_t[i] + fRand->Gaus( 0, 0.3 );
                nhith++;
            }
            
            if( Virtual_det[i] == 0 && Virtual_E[i] > 3000 && (Virtual_pdg[i] == 11 || Virtual_pdg[i] == -11 || Virtual_pdg[i] == 22) ) {
                
                p4e.SetPxPyPzE( 0.001*Virtual_px[i], 0.001*Virtual_py[i], 0.001*Virtual_pz[i], (0.001*Virtual_E[i]));
                
                p4e.SetTheta( p4e.Theta() + fRand->Gaus(0, (0.08/57.3) ) );
                
                Qp4 = kp4 - p4e;
                Rp4 = Tp4 + Qp4;
                
                double energy = Virtual_E[i] + fRand->Gaus( 0, 0.026*Virtual_E[i] );
                
                double phidiff = (p4e.Phi()*57.3 + fRand->Gaus(0, 0.08)) + 180 - ((p4h.Phi() *57.3 + fRand->Gaus(0, 0.15)));
                if( phidiff > 300 )
                    phidiff -= 360.;
                
                Double_t the     = p4e.Theta()+ fRand->Gaus(0, (0.08/57.3));
                Double_t pexp_th = 6600 / (1.0 + 6600/938*(1.0-cos(the)));
                Double_t p       = energy;
                Double_t pdiff   = (p - pexp_th)/1000.;
                
                timee = Virtual_t[i] + fRand->Gaus( 0, 0.3 );
                
                double thdiff = 57.3*(Rp4.Theta() - p4h.Theta() + fRand->Gaus(0, (0.15/57.3)));
                
                if( Event_flag == 0 ) {
                    hEarm_Ee->Fill( energy, Event_weight );
                    hEarm_E->Fill( energy, Event_weight );
                }
                else if( Event_flag == 1 ) {
                    hEarm_Eq->Fill( energy, Event_weight );
                }
                else if( Event_flag == 2 ) {
                    hEarm_Ei->Fill( energy, Event_weight );
                    hEarm_E->Fill( energy, Event_weight );
                }
                else if( Event_flag == 3 ) {
                    hEarm_Ep->Fill( energy, Event_weight );
                }
                else if( Event_flag == 4 ) {
                    hEarm_Ei->Fill( energy, Event_weight );
                    hEarm_E->Fill( energy, Event_weight );
                }
                
                
                if( goodVh && Virtual_E[i] > 4500. ) {
                    
                    double tdiff = thdiff;
                    
                    hEarmC_PhDiff->Fill( phidiff );
                    hEarmC_tDiff->Fill( tdiff );
                    
                    hEarmC4_E->Fill( energy, Event_weight );
                    
                    if( Event_flag == 0 ) {
                        hEarmC_E->Fill( energy, Event_weight );
                        hEarmC_Ee->Fill( energy, Event_weight );
                        hEarmC4_Ee->Fill( energy, Event_weight );
                        hEarmC_PhDiffe->Fill( phidiff );
                        hEarmC_tDiffe->Fill( tdiff );
                        hEarmC4_Q2->Fill( -Qp4.M2(), Event_weight );
                        hEarmC4_Th->Fill( 57.3*the, Event_weight );
                    }
                    else if( Event_flag == 1 ) {
                        hEarmC_Eq->Fill( energy, Event_weight );
                        hEarmC4_Eq->Fill( energy, Event_weight );
                        hEarmC_PhDiffq->Fill( phidiff );
                        hEarmC_tDiffq->Fill( tdiff );
                    }
                    else if( Event_flag == 2 ) {
                        hEarmC_E->Fill( energy, Event_weight );
                        hEarmC_Ei->Fill( energy, Event_weight );
                        hEarmC4_Ei->Fill( energy, Event_weight );
                        hEarmC_PhDiffi->Fill( phidiff );
                        hEarmC_tDiffi->Fill( tdiff );
                    }
                    else if( Event_flag == 3 ) {
                        hEarmC_Ep->Fill( energy, Event_weight );
                        hEarmC4_Ep->Fill( energy, Event_weight );
                        hEarmC_PhDiffp->Fill( phidiff );
                        hEarmC_tDiffp->Fill( tdiff );
                    }
                    else if( Event_flag == 4 ) {
                        hEarmC_E->Fill( energy, Event_weight );
                        hEarmC_Ei->Fill( energy, Event_weight );
                        hEarmC4_Ei->Fill( energy, Event_weight );
                        hEarmC_PhDiffi->Fill( phidiff );
                        hEarmC_tDiffi->Fill( tdiff );
                    }
                    
                    goodVe = true;
                    
                    
                    if( Virtual_E[i] > 5000. && fabs(phidiff)<0.6 && fabs(pdiff)<0.6 && fabs( (the*57.3) - 15.5) < 0.58 ) {
                        
                        
                        hEarmC5_E->Fill( energy, Event_weight );
                        
                        if( Event_flag == 0 ) {
                            hEarmC5_Ee->Fill( energy, Event_weight );
                            hEarmC5_Q2->Fill( -Qp4.M2(), Event_weight );
                            hEarmC5_Th->Fill( 57.3*the, Event_weight );
                        }
                        else if( Event_flag == 1 ) {
                            hEarmC5_Eq->Fill( energy, Event_weight );
                        }
                        else if( Event_flag == 2 ) {
                            hEarmC5_Ei->Fill( energy, Event_weight );
                        }
                        else if( Event_flag == 3 ) {
                            hEarmC5_Ep->Fill( energy, Event_weight );
                        }
                        else if( Event_flag == 4 ) {
                            hEarmC5_Ei->Fill( energy, Event_weight );
                        }
                        
                    }
                }
                
            }
            
        }
        
        bool goodh   = false;
        bool goodhod = false;
        
        for( int i =0; i < Real_Nhits; i++ ) {
            
            if( Real_det[i] == 0  ) {
                
                if( ApplyThresh && Real_edep[i] > Earm_threshold ) {
                    if( ApplyWindow && Real_t[i] < Earm_window ) {
                        
                        Int_t id = ((Real_row[i]+1)*(Real_col[i]+1));
                        hRealEarm_Edet[id]->Fill( Real_edep[i], Event_weight );
                        
                        hRealEarm_N->Fill( (Double_t)id, Event_weight );
                        hRealEarm_E->Fill( Real_edep[i], Event_weight );
                        
                        if( Event_flag == 999 ) {
                            hRealEarm_Nb->Fill( (Double_t)id, Event_weight );
                            hRealEarm_Eb->Fill( Real_edep[i], Event_weight );
                        }
                        else if( Event_flag == 0 ) {
                            hRealEarm_Ne->Fill( (Double_t)id, Event_weight );
                            hRealEarm_Ee->Fill( Real_edep[i], Event_weight );
                        }
                        else if( Event_flag == 2 ) {
                            hRealEarm_Nd->Fill( (Double_t)id, Event_weight );
                            hRealEarm_Ed->Fill( Real_edep[i], Event_weight );
                        }
                        hRealEarm_xy->Fill( Real_x[i], Real_y[i], Event_weight );
                        hRealEarm_z->Fill( Real_z[i], Event_weight );
                        hRealEarm_t->Fill( (Double_t)Real_t[i], Event_weight );
                    }
                }
            }
            
            if( goodVe && Real_det[i] == 1 ) {
                if( ApplyThresh && Real_edep[i] > Harm_threshold ) {
                    if( ApplyWindow && Real_t[i] < Harm_window ) {
                        Int_t id = ((Real_row[i]+1)*(Real_col[i]+1));
                        
                        goodh = true;
                        
                        hRealHarm_Edet[id]->Fill( Real_edep[i], Event_weight );
                        
                        hRealHarm_N->Fill( (Double_t)id, Event_weight );
                        hRealHarm_E->Fill( Real_edep[i], Event_weight );
                        
                        if( Event_flag == 999 ) {
                            hRealHarm_Nb->Fill( (Double_t)id, Event_weight );
                            hRealHarm_Eb->Fill( Real_edep[i], Event_weight );
                        }
                        else if( Event_flag == 0 ) {
                            hRealHarm_Ne->Fill( (Double_t)id, Event_weight );
                            hRealHarm_Ee->Fill( Real_edep[i], Event_weight );
                        }
                        else  if( Event_flag == 2 ) {
                            hRealHarm_Nd->Fill( (Double_t)id, Event_weight );
                            hRealHarm_Ed->Fill( Real_edep[i], Event_weight );
                        }
                        hRealHarm_xy->Fill( Real_x[i], Real_y[i], Event_weight );
                        hRealHarm_z->Fill( Real_z[i], Event_weight );
                        hRealHarm_t->Fill( (Double_t)Real_t[i], Event_weight );
                    }
                }
            }
            
            if( Real_det[i] == 2  ) {
                
                
                if( ApplyThresh && Real_edep[i] > Hodo_threshold ) {
                    if( ApplyWindow && Real_t[i] < Hodo_window ) {
                        
                        goodhod = true;
                        
                        Int_t id = ((Real_row[i]+1)*(Real_col[i]+1));
                        hRealHodo_Edet[id]->Fill( Real_edep[i], Event_weight );
                        
                        hRealHodo_N->Fill( (Double_t)id, Event_weight );
                        hRealHodo_E->Fill( Real_edep[i], Event_weight );
                        
                        if( Event_flag == 999 ) {
                            hRealHodo_Nb->Fill( (Double_t)id, Event_weight );
                            hRealHodo_Eb->Fill( Real_edep[i], Event_weight );
                        }
                        else if( Event_flag == 0 ) {
                            hRealHodo_Ne->Fill( (Double_t)id, Event_weight );
                            hRealHodo_Ee->Fill( Real_edep[i], Event_weight );
                        }
                        else  if( Event_flag == 2 ) {
                            hRealHodo_Nd->Fill( (Double_t)id, Event_weight );
                            hRealHodo_Ed->Fill( Real_edep[i], Event_weight );
                        }
                        hRealHodo_xy->Fill( Real_x[i], Real_y[i], Event_weight );
                        hRealHodo_z->Fill( Real_z[i], Event_weight );
                        hRealHodo_t->Fill( (Double_t)Real_t[i], Event_weight );
                    }
                }
            }
            
        }
        
        if (goodh ) ngoodh++;
        if (goodhod ) ngoodhod++;
        
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TLatex* tex;
    
    TCanvas* cReal1 = new TCanvas("cReal1","",1200,800);
    
    Double_t det[NHodo], sumedep[NHodo], sumdose[NHodo] ;
    Double_t edet[NHodo], esumedep[NHodo], esumdose[NHodo] ;
    
    TGraphErrors* gEdep, *gDose;
    TF1* meanedeprate, *meandoserate;
    
    cReal1->Divide(3,2);
    
    cReal1->cd(1);
    
    tex = new TLatex( 0.09, 0.7, "OpNovice2");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.5, "EArm Singles Rates");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", 0.5));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    cReal1->cd(2)->SetLogy(1);
    hRealEarm_E->SetLineColor(1);
    hRealEarm_E->Draw("hist e4");
    hRealEarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");
    hRealEarm_Ee->SetLineColor(4);
    hRealEarm_Ee->Draw("hist e4 same");
    hRealEarm_Ed->SetLineColor(2);
    hRealEarm_Ed->Draw("hist e4 same");
    hRealEarm_Eb->SetLineColor(8);
    hRealEarm_Eb->Draw("hist e4 same");
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealEarm_E->GetMean()));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal1->cd(3)->SetLogy(1);
    
    hRealEarm_N->SetLineColor(1);
    hRealEarm_N->Draw("");
    hRealEarm_N->GetYaxis()->SetTitle("Hit Rate [Hz]");
    hRealEarm_N->GetXaxis()->SetTitle("Detector ID");
    hRealEarm_N->GetYaxis()->SetRangeUser(1e0, 1.3*hRealEarm_N->GetMaximum());
    hRealEarm_Ne->SetLineColor(4);
    hRealEarm_Ne->Draw("same");
    hRealEarm_Nd->SetLineColor(2);
    hRealEarm_Nd->Draw("same");
    hRealEarm_Nb->SetLineColor(8);
    hRealEarm_Nb->Draw("same");
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", hRealEarm_N->Integral(0,NEarm)/NEarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    //  tex->Draw();
    
    cReal1->cd(4)->SetLogy(1);
    
    for(Int_t i=0; i<NEarm; i++) {
        det[i]     = i;
        sumedep[i] = (3600*hRealEarm_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
        sumdose[i] = (sumedep[i] * 100 )/(2*2*20*8/1000.); // divide by detector mass in kg then convert J/s to rad/hr
        edet[i]     = 0.;
        esumedep[i] = 0.;
        
    }
    
    gEdep = new TGraphErrors( NEarm, det, sumedep, edet, esumedep );
    gEdep->SetMarkerColor( 1 );
    gEdep->SetLineColor( 1 );
    gEdep->SetMarkerStyle( 20 );
    gEdep->SetMarkerSize( 0.4 );
    gEdep->Draw("AP");
    gEdep->GetXaxis()->SetTitle( "Detector ID");
    gEdep->GetYaxis()->SetTitle( "Edep Rate [J/hr]");
    gEdep->GetXaxis()->SetRangeUser(0,NEarm);
    gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NEarm,gEdep->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/hr", TMath::Mean(NEarm,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal1->cd(5)->SetLogy(1);
    
    gDose = new TGraphErrors( NEarm, det, sumdose, edet, esumdose );
    gDose->SetMarkerColor( 1 );
    gDose->SetLineColor( 1 );
    gDose->SetMarkerStyle( 20 );
    gDose->SetMarkerSize( 0.4 );
    gDose->Draw("AP");
    gDose->GetXaxis()->SetTitle( "Detector ID");
    gDose->GetYaxis()->SetTitle( "Dose Rate [rad/hr]");
    gDose->GetXaxis()->SetRangeUser(0,NEarm);
    gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NEarm,gDose->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", TMath::Mean(NEarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal1->cd(6);
    
    tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.4, Form("Mean edep rate = %3.2e J/hr", TMath::Mean(NEarm,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Mean dose rate = %3.2e rad/hr", TMath::Mean(NEarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.2, Form("Total dose = %3.2e rad", 720*TMath::Mean(NEarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.8, Form("Total hit rate = %3.2e Hz", hRealEarm_N->Integral(0,NEarm)/NEarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.7, Form("Beam bg hit rate = %3.2e Hz", hRealEarm_Nb->Integral(0,NEarm)/NEarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(8);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.6, Form("DIS hit rate = %3.2e Hz", hRealEarm_Nd->Integral(0,NEarm)/NEarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.5, Form("Elastic hit rate = %3.2e Hz", hRealEarm_Ne->Integral(0,NEarm)/NEarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TCanvas* cReal2 = new TCanvas("cReal2","",1200,800);
    
    cReal2->Divide(3,2);
    
    cReal2->cd(1);
    
    tex = new TLatex( 0.09, 0.7, "OpNovice2");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.5, "HArm Singles Rates");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", 0.5));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    cReal2->cd(2)->SetLogy(1);
    hRealHarm_E->SetLineColor(1);
    hRealHarm_E->Draw("hist e4");
    hRealHarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");
    hRealHarm_Ee->SetLineColor(4);
    hRealHarm_Ee->Draw("hist e4 same");
    hRealHarm_Ed->SetLineColor(2);
    hRealHarm_Ed->Draw("hist e4 same");
    hRealHarm_Eb->SetLineColor(8);
    hRealHarm_Eb->Draw("hist e4 same");
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealHarm_E->GetMean()));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal2->cd(3)->SetLogy(1);
    hRealHarm_N->SetLineColor(1);
    hRealHarm_N->Draw("");
    hRealHarm_N->GetYaxis()->SetTitle("Hit Rate [Hz]");
    hRealHarm_N->GetXaxis()->SetTitle("Detector ID");
    hRealHarm_N->GetYaxis()->SetRangeUser(1e0, 1.3*hRealHarm_N->GetMaximum());
    hRealHarm_Ne->SetLineColor(4);
    hRealHarm_Ne->Draw("same");
    hRealHarm_Nd->SetLineColor(2);
    hRealHarm_Nd->Draw("same");
    hRealHarm_Nb->SetLineColor(8);
    hRealHarm_Nb->Draw("same");
    
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", hRealHarm_N->Integral(0,NHarm)/NHarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    //  tex->Draw();
    
    
    cReal2->cd(4)->SetLogy(1);
    
    for(Int_t i=0; i<NHarm; i++) {
        det[i]     = i;
        sumedep[i] = (3600*hRealHarm_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
        sumdose[i] = (sumedep[i] * 100 )/(15*15*100*4/1000.); // divide by detector mass in kg then convert J/s to rad/hr
        edet[i]     = 0.;
        esumedep[i] = 0.;
        
    }
    
    gEdep = new TGraphErrors( NHarm, det, sumedep, edet, esumedep );
    gEdep->SetMarkerColor( 1 );
    gEdep->SetLineColor( 1 );
    gEdep->SetMarkerStyle( 20 );
    gEdep->SetMarkerSize( 0.4 );
    gEdep->Draw("AP");
    gEdep->GetXaxis()->SetTitle( "Detector ID");
    gEdep->GetYaxis()->SetTitle( "Edep Rate [J/hr]");
    gEdep->GetXaxis()->SetRangeUser(0,NHarm);
    gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHarm,gEdep->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/hr", TMath::Mean(NHarm,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal2->cd(5)->SetLogy(1);
    
    gDose = new TGraphErrors( NHarm, det, sumdose, edet, esumdose );
    gDose->SetMarkerColor( 1 );
    gDose->SetLineColor( 1 );
    gDose->SetMarkerStyle( 20 );
    gDose->SetMarkerSize( 0.4 );
    gDose->Draw("AP");
    gDose->GetXaxis()->SetTitle( "Detector ID");
    gDose->GetYaxis()->SetTitle( "Dose Rate [rad/hr]");
    gDose->GetXaxis()->SetRangeUser(0,NHarm);
    gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHarm,gDose->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", TMath::Mean(NHarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal2->cd(6);
    
    tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.4, Form("Mean edep rate = %3.2e J/hr", TMath::Mean(NHarm,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Mean dose rate = %3.2e rad/hr", TMath::Mean(NHarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.2, Form("Total dose = %3.2e rad", 720*TMath::Mean(NHarm,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.8, Form("Total hit rate = %3.2e Hz", hRealHarm_N->Integral(0,NHarm)/NHarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.7, Form("Beam bg hit rate = %3.2e Hz", hRealHarm_Nb->Integral(0,NHarm)/NHarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(8);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.6, Form("DIS hit rate = %3.2e Hz", hRealHarm_Nd->Integral(0,NHarm)/NHarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.5, Form("Elastic hit rate = %3.2e Hz", hRealHarm_Ne->Integral(0,NHarm)/NHarm));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TCanvas* cReal3 = new TCanvas("cReal3","",1200,800);
    
    cReal3->Divide(3,2);
    
    cReal3->cd(1);
    
    tex = new TLatex( 0.09, 0.7, "OpNovice2");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    
    tex = new TLatex( 0.09, 0.5, "Hodo Singles Rates");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", Hodo_threshold));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    cReal3->cd(2)->SetLogy(1);
    hRealHodo_E->SetLineColor(1);
    hRealHodo_E->Draw("hist e4");
    hRealHodo_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");
    hRealHodo_Ee->SetLineColor(4);
    hRealHodo_Ee->Draw("hist e4 same");
    hRealHodo_Ed->SetLineColor(2);
    hRealHodo_Ed->Draw("hist e4 same");
    hRealHodo_Eb->SetLineColor(8);
    hRealHodo_Eb->Draw("hist e4 same");
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealHodo_E->GetMean()));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal3->cd(3)->SetLogy(1);
    hRealHodo_N->SetLineColor(4);
    hRealHodo_N->Draw("");
    hRealHodo_N->GetYaxis()->SetTitle("Hit Rate [Hz]");
    hRealHodo_N->GetXaxis()->SetTitle("Detector ID");
    hRealHodo_N->GetYaxis()->SetRangeUser(1e-3, 1.3*hRealHodo_N->GetMaximum());
    hRealHodo_Ne->SetLineColor(4);
    hRealHodo_Ne->Draw("same");
    hRealHodo_Nd->SetLineColor(2);
    hRealHodo_Nd->Draw("same");
    hRealHodo_Nb->SetLineColor(8);
    hRealHodo_Nb->Draw("same");
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", hRealHodo_N->Integral(0,NHodo)/NHodo));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    //  tex->Draw();
    
    cReal3->cd(4)->SetLogy(1);
    
    for(Int_t i=0; i<NHodo; i++) {
        det[i]     = i;
        sumedep[i] = (3600*hRealHodo_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
        sumdose[i] = (sumedep[i] * 100 )/(3*3*10*1/1000.); // divide by detector mass in kg then convert J/s to rad/hr
        edet[i]     = 0.;
        esumedep[i] = 0.;
        
    }
    
    gEdep = new TGraphErrors( NHodo, det, sumedep, edet, esumedep );
    gEdep->SetMarkerColor( 1 );
    gEdep->SetLineColor( 1 );
    gEdep->SetMarkerStyle( 20 );
    gEdep->SetMarkerSize( 0.4 );
    gEdep->Draw("AP");
    gEdep->GetXaxis()->SetTitle( "Detector ID");
    gEdep->GetYaxis()->SetTitle( "Edep Rate [J/hr]");
    gEdep->GetXaxis()->SetRangeUser(0,NHodo);
    gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHodo,gEdep->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/hr", TMath::Mean(NHodo,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal3->cd(5)->SetLogy(1);
    
    gDose = new TGraphErrors( NHodo, det, sumdose, edet, esumdose );
    gDose->SetMarkerColor( 1 );
    gDose->SetLineColor( 1 );
    gDose->SetMarkerStyle( 20 );
    gDose->SetMarkerSize( 0.4 );
    gDose->Draw("AP");
    gDose->GetXaxis()->SetTitle( "Detector ID");
    gDose->GetYaxis()->SetTitle( "Dose Rate [rad/hr]");
    gDose->GetXaxis()->SetRangeUser(0,NHodo);
    gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHodo,gDose->GetY()) );
    
    tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", TMath::Mean(NHodo,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    cReal3->cd(6);
    
    tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.4, Form("Mean edep rate = %3.2e J/hr", TMath::Mean(NHodo,gEdep->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.3, Form("Mean dose rate = %3.2e rad/hr", TMath::Mean(NHodo,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.2, Form("Total dose = %3.2e rad", 720*TMath::Mean(NHodo,gDose->GetY()) ));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.8, Form("Total hit rate = %3.2e Hz", hRealHodo_N->Integral(0,NHodo)/NHodo));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.7, Form("Beam bg hit rate = %3.2e Hz", hRealHodo_Nb->Integral(0,NHodo)/NHodo));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(8);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.6, Form("DIS hit rate = %3.2e Hz", hRealHodo_Nd->Integral(0,NHodo)/NHodo));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    tex = new TLatex( 0.09, 0.5, Form("Elastic hit rate = %3.2e Hz", hRealHodo_Ne->Integral(0,NHodo)/NHodo));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->Draw();
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    TCanvas* cVirt1 = new TCanvas("cVirt1","",1200,800);
    cVirt1->Divide(2,1);
    
    cVirt1->cd(1)->SetLogy(1);
    
    hEarm_Ei->SetLineWidth(2) ;
    hEarm_Ei->SetLineColor(2);
    hEarm_Ei->Draw("hist e4");
    hEarm_Ei->GetYaxis()->SetRangeUser( 1e-3,  100*hEarm_E->GetMaximum());
    hEarm_Ei->GetXaxis()->SetTitle("Cluster Energy [MeV]");
    hEarm_Ei->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarm_Ee->SetLineWidth(2) ;
    hEarm_Ee->SetLineColor(4);
    hEarm_Ee->Draw("hist e4 same");
    
    hEarm_E->SetLineColor(1) ;
    hEarm_E->SetLineWidth(2) ;
    hEarm_E->Draw("hist e4 same");
    
    TLine *line = new TLine(4500, 1e-3, 4500, 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    cout << 1*hEarm_E->GetMaximum() << endl;
    
    line = new TLine(5.0, 1e-3, 5.0, 10*hEarm_E->GetMaximum());
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    //  line->Draw();
    
    tex = new TLatex(0.62,0.9,"Elastic");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.62,0.85,"Inelastic");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    tex = new TLatex(4000,2e5,"Trigger");
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(3800,7e4,"threshold");
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);
    tex->Draw();
    
    cVirt1->cd(2)->SetLogy(1);
    
    hEarm_E->SetLineColor(1);
    hEarm_E->Draw("hist e4");
    hEarm_E->GetYaxis()->SetRangeUser( 1e-3,  100*hEarm_E->GetMaximum());
    hEarm_E->GetXaxis()->SetTitle("Cluster Energy [MeV]");
    hEarm_E->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarmC_Ee->SetLineWidth(2) ;
    hEarmC_Ee->SetLineColor(4);
    hEarmC_Ee->Draw("hist e4 same");
    
    hEarmC_Ei->SetLineWidth(2) ;
    hEarmC_Ei->SetLineColor(2);
    hEarmC_Ei->Draw("hist e4 same");
    
    hEarmC_E->SetLineColor(6) ;
    hEarmC_E->SetLineWidth(2) ;
    hEarmC_E->Draw("hist e4 same");
    
    line = new TLine(4500, 1e-3, 4500, 1*30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    line = new TLine(5.0, 1e-3, 5.0, 10*hEarm_E->GetMaximum());
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    
    tex = new TLatex(0.50,0.9,"ECAL only");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.50,0.85,"ECAL+HCAL");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(6);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    
    TCanvas* cVirt2 = new TCanvas("cVirt2","",1200,800);
    cVirt2->Divide(2,2);
    
    cVirt2->cd(3)->SetLogy(1);
    
    hEarmC4_Ee->SetLineWidth(2) ;
    hEarmC4_Ee->SetLineColor(4);
    hEarmC4_Ee->Draw("hist e4");
    hEarmC4_Ee->GetYaxis()->SetRangeUser( 1e-3,  100*hEarmC4_E->GetMaximum());
    hEarmC4_Ee->GetXaxis()->SetTitle("Cluster Energy [MeV]");
    hEarmC4_Ee->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarmC4_Eq->SetLineWidth(2) ;
    hEarmC4_Eq->SetLineColor(41);
    hEarmC4_Eq->Scale(4);
    hEarmC4_Eq->Draw("hist e4 same");
    
    hEarmC4_Ep->SetLineWidth(2) ;
    hEarmC4_Ep->SetLineColor(38);
    hEarmC4_Ep->Scale(30);
    hEarmC4_Ep->Draw("hist e4 same");
    
    hEarmC4_Ei->SetLineWidth(2) ;
    hEarmC4_Ei->SetLineColor(2) ;
    hEarmC4_Ei->Draw("hist e4 same");
    
    tex = new TLatex(0.66,0.85,"Online");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    line = new TLine(5000., 1e-3, 5000., 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    cVirt2->cd(4)->SetLogy(1);
    
    hEarmC5_Ee->SetLineWidth(2) ;
    hEarmC5_Ee->SetLineColor(4);
    hEarmC5_Ee->Draw("hist e4");
    hEarmC5_Ee->GetYaxis()->SetRangeUser( 1e-3,  100*hEarmC5_E->GetMaximum());
    hEarmC5_Ee->GetXaxis()->SetTitle("Cluster Energy [MeV]");
    hEarmC5_Ee->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarmC5_Eq->SetLineWidth(2) ;
    hEarmC5_Eq->SetLineColor(41);
    hEarmC5_Eq->Scale(4);
    hEarmC5_Eq->Draw("hist e4 same");
    
    hEarmC5_Ep->SetLineWidth(2) ;
    hEarmC5_Ep->SetLineColor(38);
    hEarmC5_Ep->Scale(30);
    hEarmC5_Ep->Draw("hist e4 same");
    
    hEarmC5_Ei->SetLineWidth(2) ;
    hEarmC5_Ei->SetLineColor(2) ;
    hEarmC5_Ei->Draw("hist e4 same");
    
    tex = new TLatex(0.66,0.85,"Offline");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    cVirt2->cd(1)->SetLogy(1);
    
    hEarmC_PhDiffe->SetLineWidth(2) ;
    hEarmC_PhDiffe->SetLineColor(4);
    hEarmC_PhDiffe->Draw("hist e4");
    hEarmC_PhDiffe->GetYaxis()->SetRangeUser( 1e-3,  100*hEarmC_PhDiff->GetMaximum());
    hEarmC_PhDiffe->GetXaxis()->SetTitle("#phi_{ECAL} - #phi_{Hscint} [degrees]");
    hEarmC_PhDiffe->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarmC_PhDiffq->SetLineWidth(2);
    hEarmC_PhDiffq->SetLineColor(41);
    hEarmC_PhDiffq->Scale(4);
    hEarmC_PhDiffq->Draw("hist e4 same");
    
    hEarmC_PhDiffp->SetLineWidth(2);
    hEarmC_PhDiffp->SetLineColor(2);
    hEarmC_PhDiffp->Scale(30);
    hEarmC_PhDiffp->Draw("hist e4 same");
    
    
    hEarmC_PhDiffi->SetLineWidth(2);
    hEarmC_PhDiffi->SetLineColor(38) ;
    hEarmC_PhDiffi->Draw("hist e4 same");
    
    line = new TLine(-0.6, 1e-3, -0.6, 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    line = new TLine(0.6, 1e-3, 0.6, 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    tex = new TLatex(0.56,0.9,"Elastic");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.56,0.85,"Inelastic");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.56,0.8,"Photoproduction");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(38);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.56,0.75,"Quasi-elastic");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(41);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    cVirt2->cd(2)->SetLogy(1);
    
    hEarmC_tDiffe->SetLineWidth(2) ;
    hEarmC_tDiffe->SetLineColor(4);
    hEarmC_tDiffe->Draw("hist e4");
    hEarmC_tDiffe->GetYaxis()->SetRangeUser( 1e-3,  100*hEarmC_tDiff->GetMaximum());
    hEarmC_tDiffe->GetXaxis()->SetTitle("#theta_{Hscint} - #theta_{predicted} [degrees]");
    hEarmC_tDiffe->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    hEarmC_tDiffq->SetLineWidth(2) ;
    hEarmC_tDiffq->Scale(4);
    hEarmC_tDiffq->SetLineColor(41);
    hEarmC_tDiffq->Draw("hist e4 same");
    
    hEarmC_tDiffp->SetLineWidth(2) ;
    hEarmC_tDiffp->SetLineColor(2);
    hEarmC_tDiffp->Scale(30);
    hEarmC_tDiffp->Draw("hist e4 same");
    
    
    hEarmC_tDiffi->SetLineWidth(2) ;
    hEarmC_tDiffi->SetLineColor(38) ;
    hEarmC_tDiffi->Draw("hist e4 same");
    
    line = new TLine(-0.6, 1e-3, -0.6, 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    line = new TLine(0.6, 1e-3, 0.6, 30000);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();
    
    //-----------------------------------------------------------------------------------------------------------------------------
    
    
    TCanvas* cVirt3 = new TCanvas("cVirt3","",1200,800);
    cVirt3->Divide(2,2);
    
    cVirt3->cd(1);
    
    hEarmC4_Q2->SetLineWidth(2) ;
    hEarmC4_Q2->SetLineColor(4);
    hEarmC4_Q2->Draw("hist e4");
    hEarmC4_Q2->GetXaxis()->SetTitle("Q^2 [(GeV/c)^2]");
    hEarmC4_Q2->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    tex = new TLatex(0.58,0.85,"Online");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.58,0.78,Form("Mean = %3.2f", hEarmC4_Q2->GetMean() ));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    cVirt3->cd(2);
    
    hEarmC5_Q2->SetLineWidth(2) ;
    hEarmC5_Q2->SetLineColor(4);
    hEarmC5_Q2->Draw("hist e4");
    hEarmC5_Q2->GetXaxis()->SetTitle("Q^2 [(GeV/c)^2]");
    hEarmC5_Q2->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    tex = new TLatex(0.58,0.85,"Offline");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.58,0.78,Form("Mean = %3.2f", hEarmC5_Q2->GetMean() ));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    
    cVirt3->cd(3);
    
    hEarmC4_Th->SetLineWidth(2) ;
    hEarmC4_Th->SetLineColor(4);
    hEarmC4_Th->Draw("hist e4");
    hEarmC4_Th->GetXaxis()->SetTitle("#theta [degrees]");
    hEarmC4_Th->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    tex = new TLatex(0.58,0.85,"Online");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.58,0.78,Form("Mean = %3.2f", hEarmC4_Th->GetMean() ));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    cVirt3->cd(4);
    
    hEarmC5_Th->SetLineWidth(2) ;
    hEarmC5_Th->SetLineColor(4);
    hEarmC5_Th->Draw("hist e4");
    hEarmC5_Th->GetXaxis()->SetTitle("#theta [degrees]");
    hEarmC5_Th->GetYaxis()->SetTitle("Rate/bin [Hz]");
    
    tex = new TLatex(0.58,0.85,"Offline");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    tex = new TLatex(0.58,0.78,Form("Mean = %3.2f", hEarmC5_Th->GetMean() ));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.055);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TCanvas* cReal11 = new TCanvas("cReal11","",1200,800);
    
    hRealHarm_E->SetLineColor(1);
    hRealHarm_E->Draw("hist e4");
    hRealHarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");
    hRealHarm_Ee->SetLineColor(4);
    hRealHarm_Ee->Draw("hist e4 same");
    hRealHarm_Ed->SetLineColor(2);
    hRealHarm_Ed->Draw("hist e4 same");
    hRealHarm_Eb->SetLineColor(8);
    hRealHarm_Eb->Draw("hist e4 same");
    
}


//-----------------------------------------------------------------------------------------------------------------------------


