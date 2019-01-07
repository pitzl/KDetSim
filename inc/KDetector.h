
#ifndef _KDetector
#define _KDetector

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KDetector                                                            //
//                                                                      //
// Class for description of silicon microstrip detector                 //
// This function calculates BEGIN_LATEX                                 //
// F(x_{#frac{1}{2}}) = #prod(x < x_{#frac{1}{2}}) = #frac{1}{2}        //
// END_LATEX                                                            //
//////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "fizika.h"
#include "nrutil.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TArrayD.h"
#include "TString.h"
#include <string.h>
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TRandom3.h" // DP
#include "TMath.h"
#include "TLine.h"
#include "TF3.h"
#include "KStruct.h"
#include "TH3S.h"
#include "KGeometry.h"
#include "KMaterial.h"
#include "KField.h"
#include "TVector3.h"

class KDetector : public KGeometry, public KMaterial {

private:
  Double_t Deps;
  TRandom *ran;               //random number generator
  Double_t CalErr;            //Error of the solver
  Int_t MaxIter;              //Maximum number of iterations in eq solver
  Short_t Debug;              //Print information of drift calculation etc.

public:
  Float_t Voltage;  //Voltage
  Float_t Voltage2; //Voltage2
  TArrayF Voltages; //Array of voltages

  // Definition of space charge
  TF3  * NeffF;     //effective dopping concentration function
  TH3F * NeffH;     //effective dopping concentration histogram

  // Weigthing, electric and magnetic field
  KField Ramo[99];   // Ramo weighting potential per readout node (DP)
  KField Real;       // electric field
  Float_t B[3];      // magnetic field

  // Trapping and variables used for multiplication studies

  bool Landau;       // Landau fluct on/off
  Float_t taue;      // effective trapping time constants - used if Multiplication is ON
  Float_t tauh;      // effective trapping time constants - used if Multiplication is ON

  Int_t BreakDown;     // if break down occurs it goes to 1 otherwise is 0
  Float_t MTresh;      // treshold for taking multiplication into account
  Float_t BDTresh;     // hole multiplication - break down treshold

  // Drift parameters

  Float_t enp[3];      //entry point for the charge drift
  Float_t exp[3];      //exit point for the cahrge drift
  Int_t diff;          // Diffusion simulation (yes=1, no=0)
  Int_t average;       // Average (over how many events)
  Float_t SStep;       // Simulation step size;

  // Output histograms
  TH1F *pos;           // contribution of the holes to the total drift current
  TH1F *neg;           // contribution of the electrons  to the total drift current
  TH1F *sum;	       // total drift current
  Double_t qnode[99];  // induced charge per node (DP)

  // Constructors and destructor
  KDetector();
  ~KDetector();

  //Configuration functions
  void ResetRnd(Int_t seed) {delete ran; ran=new TRandom(seed);}; // reset the random generator
  void SetDriftHisto(Float_t x,Int_t=200);
  void SetCalculationParameters(Double_t x,Int_t y){CalErr=x; MaxIter=y;}

  // Solving the differntial equations
  void Declaration(Int_t);                 // declaration of boundary conditions
  Double_t kappa(int ,int , int , int);    // defining space charge
  Double_t V(int ,int);                    // defining voltage
  void CalField(Int_t);                    // start declaration followed by solving Poisson's equation.
  inline void CalPhyField(){CalField(0);}
  inline void CalRamoField(){CalField(1);} // should not be used for multi-electrodes (DP)

  // Calculation in case of any changes
  void SetVoltage(Float_t x,Int_t calnow=1) {Voltage=x; if(calnow) CalPhyField(); };
  void SetNeff(TF3 *neff,Int_t calnow=1) {NeffF=neff; if(calnow) CalPhyField(); };
  void SetNeff(TH3F *neff,Int_t calnow=1) {NeffH=neff; if(calnow) CalPhyField(); };

  // Simulation of drift

  void SetEntryPoint(Float_t x, Float_t y, Float_t z) {enp[0]=x; enp[1]=y; enp[2]=z;};
  void SetExitPoint(Float_t x, Float_t y, Float_t z) {exp[0]=x; exp[1]=y; exp[2]=z;};
  void MipIR(Int_t=20);
  void ShowMipIR(Int_t, Int_t=14, Int_t=1);
  void SetAverage( Int_t av ) { average = av; }; // DP

  void Drift(Double_t, Double_t, Double_t, Float_t, KStruct *, Double_t = 0);
  void CalM(KStruct *seg, Double_t *data, Int_t=-1); //multiplication calculation

  // visualization
  TH2F *Draw( std::string, Float_t=1 );
  TH1F *Draw1D( std::string, Float_t ,Int_t ,Float_t );

  // Save,read and debug
  void  Save( std::string, std::string );
  TFile *Read( std::string, std::string );
  void  SetDebug( Short_t x ) { Debug = x; };

  // precision of drift
  inline void  SetPrecision(Double_t x) { Deps = x; };
  inline Double_t GetPrecision() { return Deps; };

  ClassDef(KDetector,1)
};

#endif
