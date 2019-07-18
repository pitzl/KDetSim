#ifndef _K3D
#define _K3D

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// K3D                                                                  //
//                                                                      //
// Description of the 3D detector.                                      //
// The geometry of the                                                  //
// is defined by                                                        //
// dasdasd                                                              //
//////////////////////////////////////////////////////////////////////////


#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include "KStruct.h"
#include "TMinuit.h"
#include "KDetector.h"

class K3D : public KDetector {

private:
public:
  Int_t Col;
  Float_t CellZ;
  Float_t CellX;
  Float_t CellY;
  
  Float_t *PosD; //[Col]
  Float_t *PosX; //[Col]
  Float_t *PosY; //[Col]
  Float_t *PosR; //[Col]
  Short_t *PosW; //[Col]
  Short_t *PosM; //[Col]



  K3D(Int_t, Float_t=100, Float_t=100, Float_t=105);
  void SetUpColumn(Int_t, Float_t, Float_t, Float_t, Float_t, Short_t, Short_t);
  void SetColumnW( Int_t n, Short_t w );
  void SetUpVolume(Float_t, Float_t );
  void SetUpElectrodes(Int_t=0);
  void SetUpMaterial(Int_t mat);
  ~K3D();

  ClassDef(K3D,1) 
};


#endif

