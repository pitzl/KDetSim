{
 TF3 *f2=new TF3("f2","[0]",0,3000,0,3000,0,3000);
 f2->SetParameter(0,-2);

 KStrip *det=new KStrip(80,10,2,3,300);
 det->Voltage=400;
 det->SetUpVolume(2);

 det->SetUpElectrodes();
 det->SetBoundaryConditions();
 det->SetUpMaterial(0); //Silicon
 //det->SetUpMaterial(10); //diamond

 det->NeffF=f2;
 // det.Mat=10;
 det.SetDebug(0);
 det.CalField(0);
 det.CalField(1);

   
 det->SetEntryPoint(100,1,0.5);
 det->SetExitPoint(110,300,0.5);

 det.diff=1; 
 det->ShowMipIR(20);
}
