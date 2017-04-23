{
  Int_t j,i;
 TF3 *f3=new TF3("f2","x*y*z*0+[0]",0,3000,0,3000,0,3000);
 f2->SetParameter(0,-2);

 KStrip *det=new KStrip(240,120,1,1,150);
 det->Voltage=-200;
 det->SetUpVolume(1);
 det->SetUpElectrodes();

 //Single strip sorunded by oxide
     for(j=1;j<=det->ny;j++)
       for(i=1;i<=det->nx;i++)
	      if((i<50 || i>190) && j>130 ) det->DM->SetBinContent(i,j,1,2); else det->DM->SetBinContent(i,j,1,1);
 //
 det->SetBoundaryConditions();

 det->NeffF=f3;

 det.CalField(0);
 // det.CalField(1);
 /* det.enp[0]=60; */
/*  det.enp[1]=1; */
/*  det.exp[0]=60; */
/*  det.exp[1]=300; */
/*  det.exp[2]=0; */
/*  det.enp[2]=0; */

}
