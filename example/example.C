#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "../../missingPt/ntupler/ppTrack.C"

void example(){
 TH1D::SetDefaultSumw2(); 
 bool doPPtracking=true;
 
 TString tracking;
 if(doPPtracking) tracking="pp";
 else tracking="HI";
 TString algo=Form("ak3Calo_2014_10_21_%stracking",tracking.Data());
 
 TString infname[npthat];
 ppTrack * ftrk[npthat];
 HiTree * fhi[npthat];
 skimTree_pp * fskim[npthat];
 t * fjet[npthat];
 
 int npthat;

 int pthat[]={30,50,80,120,170,220,280,370,9999};
 double wpthat_pptracking[]={0.786183,0.0755865,0.00837536,0.00101997,0.00017521};
 double wpthat_HItracking[]={0.372166,0.0243365,0.00304998,0.000426838,4.92824e-05,8.79673e-06,1.7353e-06,2.92439e-07};

 TString directory,dataset;
 
 if(doPPtracking){ 
  npthat=5;
  for(int ipthat=0;ipthat<npthat;ipthat++){
   wpthat[i]=wpthat_pp[i];
  }
  pthat[npthat]=9999;
  directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/PYTHIA_PPForest_HIReco/";
  for(int ipthat=0;ipthat<npthat;ipthat++){
   infname[ipthat]=Form("PYTHIA_PPForest_HIReco_pthat%d_merged_forest_0",pthat[ipthat]);
  }
  dataset ="PYTHIA_PPForest_HIReco";
 }else{
  npthat=8;
  for(int ipthat=0;ipthat<npthat;ipthat++){
   wpthat[i]=wpthat_HItracking[i];
  }
  pthat[npthat]=9999;
  
  directory="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/Prod_PYTHIA_localdb_ppJEC_merged/";
  for(int ipthat=0;ipthat<npthat;ipthat++){
   infname[ipthat]=Form("HiForest_Private_PYTHIA_localdb_ppJEC_pthat%d_merged_forest_0",pthat[ipthat]);
  }
  dataset ="HiForest_Private_PYTHIA_localdb_ppJEC";
 }

 //pt bins for track efficiency correction
 int npt_eff=4;
 double ptmin_eff[]={0.5, 1, 3,  8};
 double ptmax_eff[]={  1, 3, 8,300};
 
 int npt_fake=4;
 double ptmin_fake[]={0.5, 1, 3,  8};
 double ptmax_fake[]={  1, 3, 8,300}; 
 TFile *f_multrec;
 f_multrec= new TFile(Form("../%s/multRec/multiplereco_pp.root",algo.Data()));
 TFile *f_secondary;
 f_secondary= new TFile(Form("../%s/secondary/secondary_pp.root",algo.Data()));

 TH2D * hsecondary = (TH2D*)f_secondary->Get("hpt_eta"); 
 TH2D * hmultrec = (TH2D*)f_multrec->Get("hpt_eta");

 //getting histograms for track efficiency correction 
 TFile *f_eff[npt_eff];
 TProfile2D *p_eff_accept[npt_eff];  
 TProfile *p_eff_pt[npt_eff]; 
 TProfile *p_eff_rmin[npt_eff]; 
 for(int ipt=0; ipt<npt_eff;ipt++){
   f_eff[ipt]= new TFile(Form("../%s/eff/eff_pt%d_%d_%s_dogenjet0.root",algo.Data(),(int)ptmin_eff[ipt],(int)ptmax_eff[ipt],algo.Data()));
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 TFile *f_fake[npt_fake];
 TProfile2D *p_fake_accept[npt_fake]; 
 TProfile *p_fake_pt[npt_fake]; 
 TProfile *p_fake_rmin[npt_fake]; 
 for(int ipt=0; ipt<npt_fake;ipt++){
   f_fake[ipt]= new TFile(Form("../%s/fake/fake_pt%d_%d_%s_dogenjet0.root",algo.Data(),(int)ptmin_eff[ipt],(int)ptmax_eff[ipt],algo.Data()));   
   p_fake_pt[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_pt");
   p_fake_accept[ipt]=(TProfile2D*)f_fake[ipt]->Get("p_fake_acceptance");
   p_fake_rmin[ipt]=(TProfile*)f_fake[ipt]->Get("p_fake_rmin");
 }

 TFile *outf= new TFile(Form("track_ntuple_%s_%s_testforfake_20140319.root",infname[0].Data(),algo.Data()),"recreate");
 std::string particleVars="pt:matchedpt:eta:phi:rmin:trackselect:eff:pNRec:weight";

 TNtuple *nt_particle = new TNtuple("nt_particle","",particleVars.data());
 std::string trackVars="pt:eta:phi:rmin:trackselect:trackstatus:eff:trkfake:fake:weight:multrec:secondary";

 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 //loop over events
 for(int ifile=0; ifile<npthat;ifile++){

 ftrk[ifile] = new ppTrack(Form("%s/%s.root",directory.Data(),infname[ifile].Data()));
 fhi[ifile] = new HiTree(Form("%s/%s.root",directory.Data(),infname[ifile].Data()));
 fskim[ifile] = new skimTree_pp(Form("%s/%s.root",directory.Data(),infname[ifile].Data()));
 fjet[ifile] = new t(Form("%s/%s.root",directory.Data(),infname[ifile].Data()),algo);
 int nentries = ftrk[ifile]->GetEntriesFast();
 for(int jentry=0;jentry<nentries;jentry++){
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk[ifile]->GetEntry(jentry);
  fhi[ifile]->GetEntry(jentry);
  fskim[ifile]->GetEntry(jentry);
  fjet[ifile]->GetEntry(jentry);
  if(!(fskim[ifile]->pPAcollisionEventSelectionPA && fabs(fhi[ifile]->vz)<15)) continue;
  //loop over tracks
  float weight=0;

  for(int ipthat=0;ipthat<npthat;ipthat++){
   if(fjet[ifile]->pthat<pthat[ipthat+1] && fjet[ifile]->pthat>=pthat[ipthat])weight=wpthat[ipthat];
  }
  for(int itrk=0;itrk<ftrk[ifile]->nParticle;itrk++){
    float trackselect=(ftrk[ifile]->mtrkQual[itrk] && fabs(ftrk[ifile]->mtrkDxy1[itrk]/ftrk[ifile]->mtrkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->mtrkDz1[itrk]/ftrk[ifile]->mtrkDzError1[itrk])<3 && (ftrk[ifile]->mtrkPtError[itrk]/ftrk[ifile]->mtrkPt[itrk])<0.1);
   float eta=ftrk[ifile]->pEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   float pt=ftrk[ifile]->pPt[itrk];
   if(pt<0.5) continue;
   float mpt=ftrk[ifile]->mtrkPt[itrk];
   float phi=ftrk[ifile]->pPhi[itrk];
   float rmin=99;

   
   for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
    if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
    float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }

   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_rmin=1;

   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<3)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 

   float eff=eff_accept*eff_pt*eff_rmin;
   
   //fill in the output tree
   float entry[]={pt,mpt,eta,phi,rmin,trackselect,eff,ftrk[ifile]->pNRec[itrk],weight};
   nt_particle->Fill(entry);
  }
  
  for(int itrk=0;itrk<ftrk[ifile]->nTrk;itrk++){
   float trackselect=(ftrk[ifile]->highPurity[itrk] && fabs(ftrk[ifile]->trkDxy1[itrk]/ftrk[ifile]->trkDxyError1[itrk])<3.0 && fabs(ftrk[ifile]->trkDz1[itrk]/ftrk[ifile]->trkDzError1[itrk])<3 && (ftrk[ifile]->trkPtError[itrk]/ftrk[ifile]->trkPt[itrk])<0.1);
   // float trackselect=(ftrk[ifile]->highPurity[itrk] && (ftrk[ifile]->trkPtError[itrk]/ftrk[ifile]->trkPt[itrk])<0.1);
   float eta=ftrk[ifile]->trkEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker   
   float pt=ftrk[ifile]->trkPt[itrk];

   if(pt<0.5) continue; //acceptance of the tracker
 //acceptance of the tracker
   float phi=ftrk[ifile]->trkPhi[itrk];
   float trkfake=ftrk[ifile]->trkFake[itrk];
   float trackstatus=ftrk[ifile]->trkStatus[itrk];
   float rmin=99;
 
   for(int ijet=0;ijet<fjet[ifile]->nref;ijet++){
    if(fabs(fjet[ifile]->jteta[ijet])>2 || fjet[ifile]->jtpt[ijet]<50) continue;
    float r_reco=sqrt(pow(eta-fjet[ifile]->jteta[ijet],2)+pow(acos(cos(phi-fjet[ifile]->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   float eff_accept=1;
   float eff_pt=1;
   float eff_rmin=1;
   
   float fake_pt,fake_accept,fake_rmin;
   fake_pt=fake_accept=fake_rmin=0;
   
   for(int ipt=0;ipt<npt_eff;ipt++){
    if(pt>=ptmin_eff[ipt] && pt<ptmax_eff[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));
     }     
   } 
   
   for(int ipt=0;ipt<npt_fake;ipt++){
    if(pt>=ptmin_fake[ipt] && pt<ptmax_fake[ipt]){
      fake_pt=p_fake_pt[ipt]->GetBinContent(p_fake_pt[ipt]->FindBin(pt));
      fake_accept=p_fake_accept[ipt]->GetBinContent(p_fake_accept[ipt]->GetXaxis()->FindBin(phi),p_fake_accept[ipt]->GetYaxis()->FindBin(eta));
      fake_rmin=p_fake_rmin[ipt]->GetBinContent(p_fake_rmin[ipt]->FindBin(rmin));
     }     
   }
   
   float eff=1;
   eff=eff_accept*eff_pt*eff_rmin; 
   float fake=0;
   if(pt<300)fake=fake_accept+fake_pt+fake_rmin;
   if(fake<0) fake=0;
   
   if(eff==0){
    cout<<"zero efficiency"<<" eta="<<eta<<" pt="<<pt<<" phi="<<phi<<endl;
	  if(pt>100)eff=0.8;
	  else eff=1;
   }
   
   float multrec=0;
   multrec=hmultrec->GetBinContent(hmultrec->FindBin(pt,eta));
   float secondary=0;
   secondary=hsecondary->GetBinContent(hsecondary->FindBin(pt,eta));
   
   float correction=(1-secondary)*(1-fake)/((eff)*(1+multrec));
   //fill in the output tree
   float entry[]={pt,eta,phi,rmin,trackselect,trackstatus,eff,trkfake,fake,weight,multrec,secondary};
   nt_track->Fill(entry);
  }
 }
 }
 outf->cd();
 nt_track->Write();
 nt_particle->Write();
 outf->Close();

}
