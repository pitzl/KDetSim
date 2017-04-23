
//  void KPad::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
//  {  
//    Neff->SetParameters(par);
//    GetField();
//    for(int i=0;i<PhyField.GetSize();i++) if(PhyField[i]>0) {f=1e5; return;};  
  
//    sum->Reset();
//    LaserV(entryp,exitp,20); 
//    //el->preamp(d->sum);
//    sum->GetXaxis()->Set(200,0,100);
//    f=HDif(sum,neg,0,35);
//    sum->GetXaxis()->Set(200,0,100e-9);
//    printf("DIff2=%f\n",f);
//    return;
//  }
//  void KPad::Minimize(Int_t npar, Double_t *vstart, Double_t *step)
//  {// The z values	
 
//     TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit with a maximum of 5 params
//     gMinuit->SetFCN((void *)fcn);   
//     Double_t arglist[10];   Int_t ierflg = 0;
//     arglist[0] = 1;   
//     gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
//  // Set starting values and step sizes for parameters
//     //   static Double_t vstart[3] = {250, -1. , 1};
//     //   static Double_t step[3] = {4 , 0.2 , 0.2};
//     Char_t parname[3]; parname[0]='a'; parname[1]='x';
//     for(Int_t i=0;i<npar;i++)
//       {  
//       parname[1]=48+i;
//       gMinuit->mnparm(0, parname, vstart[i], step[i], 0,0,ierflg);
//       }
//  // Now ready for minimization step   
//     arglist[0] = 500;   
//     arglist[1] = 1.;
//     gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);// Print results
   
//     Double_t amin,edm,errdef;   
//     Int_t nvpar,nparx,icstat;
//     gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
//     gMinuit->mnprin(3,amin);
//  }
