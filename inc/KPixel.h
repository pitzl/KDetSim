
#ifndef _KPixel
#define _KPixel

#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include "KStruct.h"
#include "TMinuit.h"
#include "KDetector.h"

class KPixel : public KDetector {

 private:
 public:
  Int_t nPix;
  Float_t CellZ;
  Float_t CellX;
  Float_t CellY;

  Float_t *PSx;//[nPix]
  Float_t *PSy;//[nPix]
  Float_t *PSWx;//[nPix]
  Float_t *PSWy;//[nPix]
  Float_t *PSd;//[nPix]
  Short_t *PSW;//[nPix]

  KPixel( Int_t, Float_t=200, Float_t=50, Float_t=125 );
  ~KPixel();
  void SetUpVolume( Float_t, Float_t, Float_t);
  void SetUpPixel( Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Short_t);
  void SetUpElectrodes( Int_t=0 );
  void SetPixelW( Int_t, Short_t ); // DP

  ClassDef(KPixel,1);

};

#endif

