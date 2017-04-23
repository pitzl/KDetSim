#ifndef _KImplant3D
#define _KImplant3D

#include "TMath.h"
#include "TF1.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KImplant3D                                                           //
//                                                                      //
// Calculation mesh generator                                           //
// @(#)root/html:$Id: KImplant3D.h 27910 2012-10-22 17:26:55Z Krambi $  //
// Author: Gregor Kramberger   18/10/12                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class KImplant3D {

private:
public:
  Double_t Dim[6]; // Dim[0]= width in X 
                   // Dim[1]= width in Y
                   // Dim[2]= width in Z
                   // Dim[3]= curvature in XY;
                   // Dim[4]= curvature in XZ;
                   // Dim[5]= curvature in YZ;
  TF1 *fConc;      //Concentration function
  //_______________________________________________________________________________
  KImplant3D(Double_t *,Double_t, Double_t);
  ~KImplant3D(){};

  static Double_t ImplEdge(Double_t *,Double_t *);
  static Double_t Distance(Double_t *,Double_t *, Double_t * =NULL);
  Double_t Conc(Double_t *, Double_t=-1e12);
  ClassDef(KImplant3D,1) 
};

#endif

