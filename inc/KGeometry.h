#ifndef _KGeometry
#define _KGeometry

#include "TH3I.h"
#include "TH2F.h"

TH2F * KHisProject( void *,Int_t , Int_t );
int GetNhs();

class KGeometry {

 private:

 public:

  TH3I * EG;  //electrode geometry
  TH3I * DM;  //detector material
  Int_t nx;   //x-divisions
  Int_t ny;   //y-divisions
  Int_t nz;   //z-divisions 

  // Constructors of the class
  KGeometry();
  ~KGeometry();
  KGeometry( TH3I * x ){ GetGrid( x, 0 ); };
  KGeometry( TH3I * x, TH3I * y ){ GetGrid( x, 0 ); GetGrid( x, 1 ); };

  // Setting up grids
  void GetGrid( TH3I *, Short_t = 0 ); 
  Int_t SetBoundaryConditions();
  inline int SetElecVolt(int i){ return( (i<<16)|32768 ); };
  Double_t GetStepSize( Int_t, Int_t );
  Double_t GetStepSize( Int_t, Float_t );

  // Mapping and projections
  //  TH2F *ProjectToGeometry(void *,Int_t,Int_t);
  //  Double_t InterpolateToGeometry();
  TH3F * MapToGeometry( Double_t *, Double_t = 1 );
  TH3F * GetGeom();
  void Reset( Int_t = 0, Int_t = 0 );

  // Electrodes generators
  void ElRectangle(Float_t *Pos, Float_t *Size, Int_t Wei, Int_t Mat);
  void ElLine(Float_t *r0,Float_t *r1, Float_t *W, Int_t Wei, Int_t Mat);
  void ElCylinder(Float_t *Pos,Float_t R, Float_t L,Int_t O, Int_t Wei, Int_t Mat); 

  Float_t GetLowEdge(Int_t);
  Float_t GetUpEdge(Int_t);

  ClassDef(KGeometry,1);

};

#endif
