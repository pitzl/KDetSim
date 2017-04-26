
#ifndef _KField
#define _KField

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ElectricField Class                                                  //
//                                                                      //
// Class for callculation of electric field in silicon                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "nrutil.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "TH3F.h"
#include "TH2F.h"
#include "TVector3.h"

Float_t KInterpolate2D( TH3F *, Float_t ,Float_t, Int_t=3, Int_t=1 );

//class KField : public KGeometry {
class KField {

 private:
  Int_t Method;   // Method to calculate the intermediate points (unused)
  Int_t dim;
 public:
  TH3F *U;
  TH3F *Ex;
  TH3F *Ey;
  TH3F *Ez;
  TH3F *E;

  KField() { U = NULL; Ex = NULL; Ey = NULL; Ez = NULL; E = NULL; };
  ~KField();
  void  CalFieldXYZ(Float_t x, Float_t y, Float_t z, Float_t *E); 
  TVector3 *CalFieldXYZ(Float_t x, Float_t y, Float_t z); 
  Float_t CalPotXYZ(Float_t x, Float_t y, Float_t z);
  static Float_t GetFieldPoint(Float_t *, Float_t *);
  Int_t CalField();
  Int_t GetDim(){return dim;};
  Double_t Mobility(Float_t E,Float_t T,Float_t Charg,Double_t Neff, Int_t which);
  Double_t Mobility(Float_t cx,Float_t cy,Float_t cz, Float_t T,Float_t Charg,Double_t Neff, Int_t which);
  Float_t DriftVelocity(Float_t E,Float_t Charg, Float_t T, Double_t Neff, Int_t which);
  Float_t DriftVelocity(Float_t cx,Float_t cy,Float_t cz, Float_t Charg, Float_t T, Double_t Neff, Int_t which);
  TH2F *Draw(Char_t *opt,Int_t=3, Int_t=1);
  // Double_t Alpha(Double_t , Short_t , Int_t=0);
  // inline static Double_t alpha(Double_t E) {return Alpha(E,-1,0);} 
  // inline static Double_t beta(Double_t E) {return Alpha(E,1,0);} 
  inline static Double_t alpha(Double_t E) {return 2.3e-1*TMath::Exp(-6.78*(2e7/(E*1e6)-1));};
  inline static Double_t beta(Double_t E) {return 1.3e-3*TMath::Exp(-13.2*(2e7/(E*1e6)-1));};
  Float_t M(Int_t, Float_t, Float_t, Float_t);
  //  Float_t M(TH1D *,Float_t);

  ClassDef(KField,1);

};

#endif
