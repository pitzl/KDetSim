// @(#)root/html:$Id: KMesh.h 27910 2012-10-22 17:26:55Z Krambi $
// Author: Gregor Kramberger   18/10/12

#ifndef _KMesh
#define _KMesh

#include "TMath.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KMesh                                                               //
//                                                                      //
// Calculation mesh generator                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class KMesh {

private:
  
public:
  Int_t N;
  Float_t Max;
  Float_t Min;
  //_______________________________________________________________________________
  KMesh(Float_t x0, Float_t x1=0){Max=x0; Min=x1; N=0;}
  ~KMesh(){};
  Int_t GetBins(Int_t Num, Float_t SS, Float_t ES, Float_t *X);
  Int_t GetBins(Int_t size,Float_t *Pos, Float_t *Step, Float_t *Bins);
  //  Double_t fdv() { return((Double_t) Neff1*1e-6*e_0*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
  //  Double_t fdv(Double_t Neff) {return((Double_t) Neff*1e-6*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
  ClassDef(KMesh,1) 
};

#endif


