// @(#)root/html:$Id: KStrip.h 27910 2012-10-22 17:26:55Z Krambi $
// Author: Gregor Kramberger   18/10/12

#ifndef _KStrip
#define _KStrip

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KStrip                                                               //
//                                                                      //
// Class for description of silicon microstrip detector                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "KDetector.h"

class KStrip : public KDetector {

private:
  
public:
  Float_t Pitch; //Strip pitch 
  Float_t Width; //Strip width
  Float_t Depth; //Strip depth
  Float_t CellX;  //Detector thickness
  Float_t CellY;  //Detector thickness
  Int_t   NoStrips;   //Number of strips
  Int_t   RamoStrip;  //Ramo strip

  void SetUpElectrodes(Int_t =-1);
  void SetUpVolume(Float_t x){SetUpVolume(x,x);};
  void SetUpMaterial(Int_t =1);
  void SetUpVolume(Float_t ,Float_t);
  void SetUpVolume(Float_t ,Float_t, Int_t);
  void SetUpVolume(Int_t ,Float_t *, Int_t , Float_t *);
  //_______________________________________________________________________________
  KStrip(Float_t=80, Float_t=20, Float_t=2, Int_t=3, Float_t=300);
  ~KStrip(){};
  //  Double_t fdv() { return((Double_t) Neff1*1e-6*e_0*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
  //  Double_t fdv(Double_t Neff) {return((Double_t) Neff*1e-6*TMath::Power(Thickness-1,2)/(2*perm*perm0));};
  ClassDef(KStrip,1) 
};

#endif


