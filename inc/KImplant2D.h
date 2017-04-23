// @(#)root/html:$Id: KImplant2D.h 27910 2012-10-22 17:26:55Z Krambi $
// Author: Gregor Kramberger   18/10/12

#ifndef _KImplant2D
#define _KImplant2D

#include "TMath.h"
#include "TF1.h"
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KImplant2D                                                               //
//                                                                      //
// Calculation mesh generator                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class KImplant2D {

private:
  
public:

  Double_t Dim[3]; //Dimmensions of the Implant ([0]=Size X, [1]=Size Y, [2]=Curv
  TF1 *fConc;      //Concentration function
  static Double_t ImplEdge(Double_t *,Double_t *);
  static inline Double_t ImplEdge(Double_t x,Double_t *par){return ImplEdge(&x,par);}
  static Double_t Distance(Double_t *,Double_t *, Double_t * =NULL);
  Double_t Distance1(Double_t *,Double_t *, Double_t * =NULL);
  static Double_t Derivative(Double_t *,Double_t *);
  static inline Double_t Derivative(Double_t x,Double_t *par){return Derivative(&x,par);}
  inline Double_t PDistance(Double_t *x,Double_t *y) {return TMath::Sqrt(TMath::Power(x[0]-y[0],2)+TMath::Power(x[1]-y[1],2));}
  Double_t Conc(Double_t *, Double_t=-1e12);
  //_______________________________________________________________________________
  KImplant2D(Double_t *,Double_t, Double_t);
  ~KImplant2D(){};
  ClassDef(KImplant2D,1) 
};

#endif

