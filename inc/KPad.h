#ifndef _KPad
#define _KPad

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KPad                                                                //
//                                                                      //
// Description of the pad detector                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "fizika.h"
#include "TObject.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include "KStruct.h"
#include "TMinuit.h"
#include "KDetector.h"

Double_t laser(Double_t *, Double_t *);

class KPad : public KDetector {

private:
//Runge Kutta method for solving the field
  void           rk4(float *,float *,int,float,float,float*); 
  Float_t        rtbis(float, float, float);
  Float_t        PoEqSolve(Float_t);
  void           Derivs(float x,float *,float *);
  TArrayF PhyPot;       //electric potential
  TArrayF PhyField;     //electric field 
public:
   TF1     *Neff;   // effective dopping concentration 
   Float_t CellY;   // thickness of the diode
   Float_t CellX;   // width of the diode

   KPad(Float_t=50,Float_t=301);
  ~KPad(); 
   void SetUpVolume(Float_t);
   void SetUpElectrodes();
 
 
   TGraph   *DrawPad(char*);
   void CalField(Int_t what) {if(what==1) GetRamoField(); else GetField();}
   void GetRamoField(TH1F *rf);
   void GetField(TH1F *rf);
   void GetField(TF1 *rf);
   void    GetField();
   void    GetRamoField();

  ClassDef(KPad,1) 
};


#endif











