#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

#include "volumes.C"

#include "TRandom3.h"

using namespace std;

TFile *f_ncdeltarad = new TFile("ncdeltarad_overlay_run1_v19.4.root");
TFile *f_ncpi0 = new TFile("ncpi0_overlay_run1_v19.4.root");

TDirectory* dir_ncdeltarad = (TDirectory*)f_ncdeltarad->Get("singlephoton");
TDirectory* dir_ncpi0 = (TDirectory*)f_ncpi0->Get("singlephoton");

TTree* tree_ncdeltarad = (TTree*)dir_ncdeltarad->Get("vertex_tree");
TTree* tree_ncpi0 = (TTree*)dir_ncpi0->Get("vertex_tree");


void random_vtx(){

  TH2D * totalXY = new TH2D ("totalXY","totalXY",400,0.,400.,400,-200.,200.);
  TH2D * XYact = new TH2D ("XYact","XYact",400,0.,400.,400,-200.,200.);
  TH2D * XY5cm = new TH2D ("XY5cm","XY5cm",400,0.,400.,400,-200.,200.);
  TH2D * XY10cm = new TH2D ("XY10cm","XY10cm",400,0.,400.,400,-200.,200.);
  TH2D * XYscb = new TH2D ("XYscb","XYscb",400,0.,400.,400,-200.,200.);

  TH2D * XYscb2cm = new TH2D ("XYscb2cm","XYscb2cm",400,0.,400.,400,-200.,200.);


  TH2D * XZ10cm = new TH2D ("XZ10cm","XZ10cm",400,0.,400.,2000,-400.,1600.);
  TH2D * XZ5cm = new TH2D ("XZ5cm","XZ5cm",400,0.,400.,2000,-400.,1600.);
  TH2D * XZact = new TH2D ("XZact","XZact",400,0.,400.,2000,-400.,1600.);
  TH2D * XZscb = new TH2D ("XZscb","XZscb",400,0.,400.,2000,-400.,1600.);

  TH2D * XZscb2cm = new TH2D ("XZscb2cm","XZscb2cm",400,0.,400.,2000,-400.,1600.);

  TRandom3 *r3=new TRandom3();
  r3->SetSeed(0);

  Double_t count_5cm = 0.;
  Double_t count_scb = 0.;
  Double_t count_act = 0.;
  Double_t count_10cm = 0.;

  Double_t count_scb2cm = 0.;

  TVector3 pt(40., 40., 50.);



  //vector<TGeoPolygon *> ply_vec;
  //  Bool_t testbool = load_scb(ply_vec);

  //cout << "sizeof ply_vec  " << ply_vec.size() << endl;

  Int_t tot = 10000000;

  for (Int_t i = 0; i < tot; i++){
    //x btw 0,400
    Double_t x = r3->Rndm();
    x*=400.;
    Double_t y = r3->Rndm();
    y*=400.;
    y-=200.;
    Double_t z = r3->Rndm();
    z*=2000.;
    z-=400.;
    
    pt.SetXYZ(x,y,z);

    if (is_contained(10.,pt)) {
      count_10cm+=1.;
      XY10cm->Fill(x,y);
      XZ10cm->Fill(x,z);
    }
    if (is_contained(5.,pt)) {
      count_5cm+=1.;
      XY5cm->Fill(x,y);
      XZ5cm->Fill(x,z);
    }
    if (is_contained_scb(0.,pt)) {
      count_scb+=1.;
      XYscb->Fill(x,y);
      XZscb->Fill(x,z);
    }
    if (is_contained(0.,pt)) {
      count_act+=1.;
      XYact->Fill(x,y);
      XZact->Fill(x,z);
    }
    if (is_contained_scb(2.,pt)) {
      count_scb2cm+=1.;
      XYscb2cm->Fill(x,y);
      XZscb2cm->Fill(x,z);
    }
  }

  auto cxy10 = new TCanvas("cxy10","cxy10");
  cxy10->cd();
  XY10cm->Draw("colz");
  auto cxz10 = new TCanvas("cxz10","cxz10");
  cxz10->cd();
  XZ10cm->Draw("colz");

  auto cxyscb = new TCanvas("cxyscb","cxyscb");
  cxyscb->cd();
  XYscb->Draw("colz");
  auto cxzscb = new TCanvas("cxzscb","cxzscb");
  cxzscb->cd();
  XZscb->Draw("colz");

  auto cxyscb2cm = new TCanvas("cxyscb2cm","cxyscb2cm");
  cxyscb2cm->cd();
  XYscb2cm->Draw("colz");
  auto cxzscb2cm = new TCanvas("cxzscb2cm","cxzscb2cm");
  cxzscb2cm->cd();
  XZscb2cm->Draw("colz");

  cout << "total points: " << tot << " , active vol: "<< count_act <<" , 5cm fid : " << count_5cm <<" , 10cm fid : " << count_10cm << " , scb : " << count_scb << " , scb2cm : " << count_scb2cm<< endl; 
  cout << "total points: " << tot << " , active vol : " <<(Double_t)count_act/tot<< " , 5cm fid : " <<(Double_t)count_5cm/tot << " , 10cm fid : " <<(Double_t)count_10cm/tot << " , scb : " << (Double_t)count_scb/tot << " , scb2cm : " << (Double_t)count_scb2cm/tot << endl;


}

void cuts(Double_t fid_cut, Double_t trk_cut, Bool_t isdelta);

Double_t cuts(Double_t fid_cut, Double_t trk_cut, Bool_t isdelta, Bool_t plots);

void cuts(){
  cuts(-1,-1, true);
}

void cuts(Double_t fid_cut, Double_t trk_cut, Bool_t isdelta){
  Double_t a = cuts(-1,-1, true, true);
  cout << "final eff : " << a << endl;
}

void sce(){
  
  TVector3 mock_vec(40., 40., 50.);
  Bool_t testbool = is_contained_scb(0.,mock_vec);

  for (int i =1 ; i< 15; i++){
    mock_vec.SetZ(i*100+50.);
    Bool_t testbool = is_contained_scb(0.,mock_vec);
  }

}

void cuts2D(){

  Double_t from_edge[10] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};

  ofstream file("deltarad.txt");
  //file << "Hello, world!" << std::endl;
  //stdout = fopen ("standard-output-file", "w");
  for (int i=0;i<10;i++){
    for (int j=0; j<10; j++){
      Bool_t isdelta = true;
      Double_t eff = cuts(from_edge[i],from_edge[j],isdelta,1);
      //stdout = fopen ("standard-output-file", "w");

      file << "isdelta? " << isdelta << " , fidcut :  "<<from_edge[i] << " , trkcut : " <<from_edge[j]<<" , eff : "<< eff <<endl; 
      
    }
  }
    
}

Double_t cuts(Double_t fid_cut, Double_t trk_cut, Bool_t isdelta, Bool_t plots){

  TTree* tree;
  if(isdelta){
  tree= (TTree*)dir_ncdeltarad->Get("vertex_tree");
  }  
  else{
  tree= (TTree*)dir_ncpi0->Get("vertex_tree");
  }

  Int_t reco_vtx_num;
  Double_t reco_vertex_x, reco_vertex_y, reco_vertex_z ;

  tree->SetBranchAddress("reco_vertex_size",&reco_vtx_num);
  tree->SetBranchAddress("reco_vertex_x",&reco_vertex_x);
  tree->SetBranchAddress("reco_vertex_y",&reco_vertex_y);
  tree->SetBranchAddress("reco_vertex_z",&reco_vertex_z);

  Int_t reco_asso_tracks, reco_asso_showers;
  tree->SetBranchAddress("reco_asso_tracks",&reco_asso_tracks);
  tree->SetBranchAddress("reco_asso_showers",&reco_asso_showers);


  vector<double> *reco_track_endx=NULL;
  vector<double> *reco_track_endy=NULL;
  vector<double> *reco_track_endz=NULL;

  tree->SetBranchAddress("reco_track_endx",&reco_track_endx);
  tree->SetBranchAddress("reco_track_endy",&reco_track_endy);
  tree->SetBranchAddress("reco_track_endz",&reco_track_endz);

  Long64_t nentries = tree->GetEntries();

  //define cutflow historagm here.
  if (fid_cut<0) fid_cut =5.;
  if (trk_cut<0) trk_cut = 2.;

  char *histname = new char[100];
  sprintf(histname,"h_cuts_vtx_%.1fcm_trkend_%.1fcm",fid_cut,trk_cut);
  TH1F *h_cuts= new TH1F("h_cuts",histname,7,0,7);

  char *binname = new char[100];
  char *filename = new char[100];

  for (Long64_t i=0;i<1000;i++) {

    tree->GetEntry(i);

    h_cuts->Fill(0);
    h_cuts->GetXaxis()->SetBinLabel(1,"nocut");
      
    if(reco_vtx_num==1){
      h_cuts->Fill(1);
      h_cuts->GetXaxis()->SetBinLabel(2,"vtx==1");
    }
    else continue;

    if(reco_asso_tracks==1){
      h_cuts->Fill(2);
      h_cuts->GetXaxis()->SetBinLabel(3,"trk==1");
    }
    else continue;

    if(reco_asso_showers==1){
      h_cuts->Fill(3);
      h_cuts->GetXaxis()->SetBinLabel(4,"shw==1");
    }
    else continue;
    
    TVector3 reco_vertex_vec(reco_vertex_x, reco_vertex_y, reco_vertex_z);
    if(is_contained(fid_cut,reco_vertex_vec)){
      h_cuts->Fill(4);
      sprintf(binname,"vtx_%1fcm",fid_cut);
      h_cuts->GetXaxis()->SetBinLabel(5,binname);
    }
    if(is_contained_scb(0.,reco_vertex_vec)){
      h_cuts->Fill(5);
      //sprintf(binname,"vtx_%1fcm",fid_cut);
      h_cuts->GetXaxis()->SetBinLabel(6,"vtx_in_scb");
    }
    //    else 
    continue;

    TVector3 reco_track_end_vec(reco_track_endx->at(0), reco_track_endy->at(0), reco_track_endz->at(0));
    if(is_contained(trk_cut,reco_track_end_vec)){
      h_cuts->Fill(6);
      sprintf(binname,"trk_end_%1fcm",trk_cut);
      h_cuts->GetXaxis()->SetBinLabel(7,binname);
    }
    else continue;

  }
  auto c = new TCanvas("c","c");
  h_cuts->Draw();
  //  gPad->SetLogy(1);

  //h_cuts->SaveAs("cutflow_ncdelta.pdf");                                                                                                                  
  //  h_cuts->SaveAs("cutflow_ncpi0.pdf");


  Int_t nocut = h_cuts->GetBinContent(1);
  Int_t tpcut = h_cuts->GetBinContent(4);

  sprintf(binname,"shw==1\ntpcut=%.3f",tpcut/(Double_t)nocut); 
  h_cuts->GetXaxis()->SetBinLabel(4,binname);

  sprintf(binname,"vtx_%.1fcm_%.3f",fid_cut,h_cuts->GetBinContent(5)/(Double_t)tpcut);
  h_cuts->GetXaxis()->SetBinLabel(5,binname);

  sprintf(binname,"trkend_%.1fcm_%.3f",trk_cut,h_cuts->GetBinContent(7)/(Double_t)tpcut);
  h_cuts->GetXaxis()->SetBinLabel(7,binname);

  if(isdelta)
  sprintf(filename,"deltarad_cutflow_vtx_%.0fcm_trkend_%.0fcm.pdf",fid_cut,trk_cut);

  else{
  sprintf(filename,"ncpi0_cutflow_vtx_%.0fcm_trkend_%.0fcm.pdf",fid_cut,trk_cut);
  }

  if (plots)  c->SaveAs(filename);

  cout << "no cut : " << nocut << endl;
  cout << "af. tp cut: " << tpcut << endl;
  cout << "ver in fid: " << h_cuts->GetBinContent(5)/(Double_t)tpcut << endl;
  cout << "trk in 2cm: " << h_cuts->GetBinContent(6)/(Double_t)tpcut << endl;

  Double_t eff = h_cuts->GetBinContent(6)/(Double_t)tpcut;

  //  delete h_cuts;

  return eff;

}
/*
bool is_contained(Double_t cut, TVector3 pt){

  //Double_t cut = 5.;

  Double_t xmin = cut;
  Double_t xmax = 256-cut;
  Double_t ymin = -117+cut;
  Double_t ymax = 117-cut;
  Double_t zmin = cut;
  Double_t zmax = 1036-cut;

  //cout << "xmin "<< xmin << "xmax "<<xmax << "ymin " << ymin <<"ymax "<<ymax  << "zmin "<<zmin <<"zmax "<<zmax <<endl;

  //  cout << "pt.X()<<endl;

  bool x_contain = (xmin < pt.X()) && (pt.X() < xmax);
  bool y_contain = (ymin < pt.Y()) && (pt.Y() < ymax);
  bool z_contain = (zmin < pt.Z()) && (pt.Z() < zmax);

  //cout << "x contain"<< x_contain << endl; 
  //cout << "y contain"<<y_contain << endl;
  //cout << "z contain"<<z_contain << endl;
  return (x_contain && y_contain && z_contain);
    // (reco_vertex_x > 5 && reco_vertex_x < 256-5 && reco_vertex_y > -117+5 && reco_vertex_y < 117-5 && reco_vertex_z > 5 && reco_vertex_z < 1036-5)

}
*/
