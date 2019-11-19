#include <TGeometry.h>
#include <params.C>
#include <TMath.h>

bool is_contained(Double_t cut, TVector3 pt){

  // these numbers from /cvmfs/uboone.opensciencegrid.org/products/ubcore/v08_00_00_23/gdml/microboonevX.gdml "TPCActive"
  Double_t xmin = cut;
  Double_t xmax = 256.35-cut;
  Double_t ymin = -116.5+cut;
  Double_t ymax = 116.5-cut;
  Double_t zmin = cut;
  Double_t zmax = 1036.8-cut;

  bool x_contain = (xmin < pt.X()) && (pt.X() < xmax);
  bool y_contain = (ymin < pt.Y()) && (pt.Y() < ymax);
  bool z_contain = (zmin < pt.Z()) && (pt.Z() < zmax);

  return (x_contain && y_contain && z_contain);

}

bool load_scb(vector<TGeoPolygon*>& zpolygons){
  TGeoManager *geom = new TGeoManager("save scb", "save scb");
  
  cout << "size of " << zpolygons.size() << endl;
  Double_t tbi = -10000.; //means "to be initialized"                                                                                                                                
  //torso
  Double_t ptX[6] = {0., tbi, YX_TOP_x2_array, YX_BOT_x2_array, tbi, 0.};
  Double_t ptY[6] = {YX_TOP_y1_array, YX_TOP_y1_array, tbi, tbi, YX_BOT_y1_array, YX_BOT_y1_array};

  TGeoPolygon *polyXY = new TGeoPolygon(6);

  for (Int_t z_idx_YX=0; z_idx_YX<10; z_idx_YX++) {
    ptX[1] = YX_TOP_x1_array[z_idx_YX+1];
    ptX[4] = YX_BOT_x1_array[z_idx_YX+1];
    ptY[2] = YX_TOP_y2_array[z_idx_YX+1];
    ptY[3] = YX_BOT_y2_array[z_idx_YX+1];
    polyXY->SetXY(ptX,ptY);
    polyXY->FinishPolygon();
  zpolygons.push_back(polyXY);
  }


  cout << "size of " <<zpolygons.size() << endl;
  return 0;


}
bool is_contained_scb(Double_t cut, TVector3 pt);

Double_t dist_scb(Double_t cut, TVector3 pt);

Double_t dist_scb(TVector3 pt){
  return dist_scb(0.,pt);
}

Double_t dist_scb(Double_t cut, TVector3 pt){
  if (!is_contained_scb(cut, pt)) return -1.;

  Int_t z_idx = (Int_t)pt.Z()/100;

  Int_t z_idx_YX = z_idx;//YX-view effective z index, it is only different for z > 10m area, where we want to appliy 9m<z<10m YX boundary, still need to keep the original z_idx bc \it's needed in ZX-view                                                                                                                                                               
  if (z_idx_YX==10) z_idx_YX-=1;

  Double_t tbi = -10000.; //means "to be initialized"                                                                                                                                
  Double_t ptX[6] = {0.+cut, tbi, YX_TOP_x2_array-cut, YX_BOT_x2_array-cut, tbi, 0.+cut};
  Double_t ptY[6] = {YX_TOP_y1_array-cut, YX_TOP_y1_array-cut, tbi, tbi, YX_BOT_y1_array+cut, YX_BOT_y1_array+cut};

  TGeoPolygon *polyXY = new TGeoPolygon(6);

  ptX[1] = YX_TOP_x1_array[z_idx_YX+1];
  ptX[4] = YX_BOT_x1_array[z_idx_YX+1];
  ptY[2] = YX_TOP_y2_array[z_idx_YX+1];
  ptY[3] = YX_BOT_y2_array[z_idx_YX+1];

  polyXY->SetXY(ptX,ptY);
  polyXY->FinishPolygon();
  Double_t testpt[2] = {pt.X(), pt.Y()};

  Bool_t XY_contain = polyXY->Contains(testpt);
  if(0<z_idx && z_idx<10) {
    Double_t up_z = pt.Z()-0.; // gonna bet, if it's middle enough to be up_z or down_z is smaller than the safefy in this z regime (1m,10m), it is safe to set up_z = z-0, down_z=1036.8-z 
    Double_t down_z = 1036.8-pt.Z();
    Double_t min_z =  TMath::Min(up_z,down_z);
    
    int safety_idx = 0;
    Double_t xy_d = polyXY->Safety(testpt, safety_idx);

    Double_t min_d = TMath::Min(xy_d, min_z);

    return min_d;
  }

  //up or down
  Double_t top_y = 116.5-pt.Y();
  Double_t bottom_y = pt.Y()+116.5;
  Double_t min_y = TMath::Min(top_y, bottom_y);

  Int_t y_idx = (pt.Y()+116.)/24;
  if (pt.Y()<-116. && pt.Y()>-116.5) y_idx = 0; //just the 0.5cm                                                                                                                   
  if(y_idx<0 || y_idx>9) return -1;


  //upstream
  if(z_idx==0){
    
    Double_t ZX_Up_z1_array = 0.;
    Double_t ZX_Up_x1_array = 120.;
    Double_t ZX_Up_z2_array = 11.;
    Double_t ZX_Up_x2_array = 256.;//upstream    

    Double_t ptX_Up[5] = {0.+cut, ZX_Up_x1_array, ZX_Up_x2_array-cut, ZX_Up_x2_array-cut, 0+cut};
    Double_t ptZ_Up[5] = {0.+cut,0.+cut,ZX_Up_z2_array, 55555., 55555.};//55555 is arbitrarly far out from the upstream.

    TGeoPolygon *polyXZ_Up = new TGeoPolygon(5);
    polyXZ_Up->SetXY(ptX_Up,ptZ_Up);
    polyXZ_Up->FinishPolygon();

    Double_t testpt_Up[2] = {pt.X(), pt.Z()};

    Int_t safety_idx_Up = 0;
    Double_t xz_d_Up = polyXZ_Up->Safety(testpt_Up, safety_idx_Up);
    Double_t min_d_Up = min(xz_d_Up, min_y);

    return min_d_Up;
  }  

  //downstream                                                                                                                                                  
  if(z_idx==10){
    
    Double_t ZX_Dw_z1_array     = 1037.;
    Double_t ZX_Dw_x1_array[11] = {0., 120.00, 115.24, 108.50, 110.67, 120.90, 126.43, 140.51, 157.15, 120.00, 120.00};
    Double_t ZX_Dw_z2_array[11] = {0., 1029.00, 1029.12, 1027.21, 1026.01, 1024.91, 1025.27, 1025.32, 1027.61, 1026.00, 1026.00};
    Double_t ZX_Dw_x2_array     = 256.;//downstream              

    Double_t ptX_Dw[5] = {0.+cut,ZX_Dw_x2_array-cut, ZX_Dw_x2_array-cut, tbi, 0.+cut};
    Double_t ptZ_Dw[5] = {55.,55.,tbi,ZX_Dw_z1_array-cut, ZX_Dw_z1_array-cut};// 55 is arbitrarily far from the downstream.

    ptX_Dw[3] = ZX_Dw_x1_array[y_idx+1];
    ptZ_Dw[2] = ZX_Dw_z2_array[y_idx+1];

    TGeoPolygon *polyXZ_Dw = new TGeoPolygon(5);
    polyXZ_Dw->SetXY(ptX_Dw,ptZ_Dw);
    polyXZ_Dw->FinishPolygon();

    Double_t testpt_Dw[2] = {pt.X(), pt.Z()};
    
    int safety_idx_Dw = 0;
    Double_t xz_d_Dw = polyXZ_Dw->Safety(testpt_Dw, safety_idx_Dw);
    Double_t min_d_Dw = min(xz_d_Dw, min_y);

    return min_d_Dw;
  }

  return -1.;

}

bool is_contained_scb(TVector3 pt){
  return is_contained_scb(0., pt);
}

bool is_contained_scb(Double_t cut, TVector3 pt){

  //TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");

  if (!is_contained(0., pt)) return 0; // is it in active volume?

  Int_t z_idx = (Int_t)pt.Z()/100;

  Int_t z_idx_YX = z_idx;//YX-view effective z index, it is only different for z > 10m area, where we want to appliy 9m<z<10m YX boundary, still need to keep the original z_idx bc it's needed in ZX-view

  if (z_idx_YX==10) z_idx_YX-=1;
  /*
  Double_t YX_TOP_y1_array     = 116.;
  Double_t YX_TOP_x1_array[11] = {0., 150.00, 132.56, 122.86, 119.46, 114.22, 110.90, 115.85, 113.48, 126.36, 144.21};
  Double_t YX_TOP_y2_array[11] = {0., 110.00, 108.14, 106.77, 105.30, 103.40, 102.18, 101.76, 102.27, 102.75, 105.10};
  Double_t YX_TOP_x2_array = 256.;
    
  Double_t YX_BOT_y1_array     = -115.;
  Double_t YX_BOT_x1_array[11] = {0., 115.71, 98.05, 92.42, 91.14, 92.25, 85.38, 78.19, 74.46, 78.86, 108.90};
  Double_t YX_BOT_y2_array[11] = {0., -101.72, -99.46, -99.51, -100.43, -99.55, -98.56, -98.00, -98.30, -99.32, -104.20};
  Double_t YX_BOT_x2_array = 256.;
  */
  Double_t tbi = -10000.; //means "to be initialized"
  
  Double_t ptX[6] = {0.+cut, tbi, YX_TOP_x2_array-cut, YX_BOT_x2_array-cut, tbi, 0.+cut};
  Double_t ptY[6] = {YX_TOP_y1_array-cut, YX_TOP_y1_array-cut, tbi, tbi, YX_BOT_y1_array+cut, YX_BOT_y1_array+cut};

  TGeoPolygon *polyXY = new TGeoPolygon(6);

  ptX[1] = YX_TOP_x1_array[z_idx_YX+1];
  ptX[4] = YX_BOT_x1_array[z_idx_YX+1];
  ptY[2] = YX_TOP_y2_array[z_idx_YX+1];
  ptY[3] = YX_BOT_y2_array[z_idx_YX+1];

  polyXY->SetXY(ptX,ptY);
  polyXY->FinishPolygon();
  Double_t testpt[2] = {pt.X(), pt.Y()};

  //cout << "is testpt ("<< pt.X()<<", "<<pt.Y()<<") contrained? "<<  polyXY->Contains(testpt)<<endl;
  //cout << "area ? " << polyXY->Area()<< endl;

  //polyXY->Draw();    
  
  Bool_t XY_contain = polyXY->Contains(testpt);

  if(0<z_idx && z_idx<10) return XY_contain;

  // if z_idx==0 or z_idx==10, they need xz view,  

  /// ZX view has Y dependence: Y sub-range from -116 to 116cm per 24cm
  Double_t ZX_Up_z1_array = 0.;
  Double_t ZX_Up_x1_array = 120.;
  Double_t ZX_Up_z2_array = 11.;
  Double_t ZX_Up_x2_array = 256.;//upstream
    
  Double_t ZX_Dw_z1_array     = 1037.;
  Double_t ZX_Dw_x1_array[11] = {0., 120.00, 115.24, 108.50, 110.67, 120.90, 126.43, 140.51, 157.15, 120.00, 120.00};
  Double_t ZX_Dw_z2_array[11] = {0., 1029.00, 1029.12, 1027.21, 1026.01, 1024.91, 1025.27, 1025.32, 1027.61, 1026.00, 1026.00};
  Double_t ZX_Dw_x2_array     = 256.;//downstream

  Int_t y_idx = (pt.Y()+116.)/24;
  
  if (pt.Y()<-116. && pt.Y()>-116.5) y_idx = 0; //just the 0.5cm 

  if(y_idx<0 || y_idx>9) return 0;

  Bool_t ZX_contain = false;

  if(z_idx==0){
    Double_t ptX_Up[5] = {0.+cut, ZX_Up_x1_array, ZX_Up_x2_array-cut, ZX_Up_x2_array-cut, 0+cut};
    Double_t ptZ_Up[5] = {0.+cut,0.+cut,ZX_Up_z2_array, 100., 100.};

    TGeoPolygon *polyXZ_Up = new TGeoPolygon(5);
    polyXZ_Up->SetXY(ptX_Up,ptZ_Up);
    polyXZ_Up->FinishPolygon();

    Double_t testpt_Up[2] = {pt.X(), pt.Z()};
    ZX_contain = polyXZ_Up->Contains(testpt_Up);
  }

  else if (z_idx==10){
    Double_t ptX_Dw[5] = {0.+cut,ZX_Dw_x2_array-cut, ZX_Dw_x2_array-cut, tbi, 0.+cut};
    Double_t ptZ_Dw[5] = {1000.,1000.,tbi,ZX_Dw_z1_array-cut, ZX_Dw_z1_array-cut};

    ptX_Dw[3] = ZX_Dw_x1_array[y_idx+1];
    ptZ_Dw[2] = ZX_Dw_z2_array[y_idx+1];

    TGeoPolygon *polyXZ_Dw = new TGeoPolygon(5);
    polyXZ_Dw->SetXY(ptX_Dw,ptZ_Dw);
    polyXZ_Dw->FinishPolygon();

    Double_t testpt_Dw[2] = {pt.X(), pt.Z()};
    ZX_contain = polyXZ_Dw->Contains(testpt_Dw);
  }


  return (XY_contain && ZX_contain);


}

