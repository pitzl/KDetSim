{
  gROOT->Reset();
  gSystem->Load( "/home/pitzl/silicon/KDetSim/lib/KDetSim.sl" );
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // DEFINITION OF VARIABLES                                                            //
  ////////////////////////////////////////////////////////////////////////////////////////

  short makeScan = 0;
  short savePlot = 0;
  
  short logicConnected = 0;	// additional logic section above oxide: floating ... 0 (default), biased ... 1
  
  double entryPointX = 110;	// entry and exit point in X (for a single injection)
  double entryPointY = 70;	// entry and exit point in Y
  
  double stripSize = 0.5*50; 	   // one half strip width in microns
  double stripPitch = 100; 	   // in microns
  double substrateThickness = 100; // substrate thickness = 100 microns
  double SOIthickness = 5;	   // thickness of the logic layer above oxide in microns
  double BOXthickness = 3;	   // thickness of the buried oxide layer
  double meshSize = 1;	           // size of the calculation mesh in microns 
                                   // this macro assumes meshsize of 1 micron
                                   // it can be smaller but, some addtional things should be changed in 
                                   // the macro which  I didn't put in order to make it simpler
  
  double voltage1 = -200;	// Voltage on the backplane
  double voltage2 = -0;	        // Voltage on the logic (nWell = GND)
  
  double oxideCharge = +300;      // concentration of the oxide charge in 1e12 cm-3
  double accumulationCharge = -0; // leave at 0
  
  // Block to name the plots
  char histoTitle[128];  
  char fileName[128];
  sprintf (histoTitle, "CCE V_{bias} = %.0lf V, logic %s%s", -voltage1, !logicConnected ? "floating" : "biased", !logicConnected ? "" : Form(" %.0lf V", voltage2));
  sprintf (fileName, Form("CCE_V_%.0lfV_logic%s%s", -voltage1, !logicConnected ? "floating" : "biased", !logicConnected ? "" : Form("_%.0lfV", -voltage2)));

  ////////////////////////////////////////////////////////////////////////////////////////  
  // START OF THE INITIALIZATION                                                        //
  ////////////////////////////////////////////////////////////////////////////////////////

  KDetector det;       // intialization of the detector - empty detector object
  det.nx=3*stripPitch; // define the width of that box , nx=integer number of nodes in x
  det.ny=substrateThickness + BOXthickness + SOIthickness; // ny=integer number of nodes in y
  det.nz=1;                                                // nz=integer number of nodes in z - if 2D then nz=1
  
  ////////////////////////////////////////////////////////////////////////////////////////    
  //INITIALIZE GEOMETRY                                                                 //
  // EG is a geometry 3D histogram where you define electrodes                          //
  // in the next step the simulator uses it to calculate also the boundary conditions   //
  ////////////////////////////////////////////////////////////////////////////////////////  

  det.EG=new TH3I("EG","EG",det.nx,0,det.nx,det.ny,0,det.ny,det.nz,0,1);
  det.EG->GetXaxis()->SetTitle("x [#mum]");
  det.EG->GetYaxis()->SetTitle("y [#mum]");
  det.EG->GetZaxis()->SetTitle("z [#mum]");

  ////////////////////////////////////////////////////////////////////////////////////////    
  //INITIALIZE MATERIAL                                                                 //
  // DM is a material histogram                                                         //
  // basically for silicon detectors: silicon, PolySi, oxide, Al, air are the options   //
  ////////////////////////////////////////////////////////////////////////////////////////  

  det.DM=new TH3I("DM","DM",det.nx,0,det.nx,det.ny,0,det.ny,det.nz,0,1);
  det.DM->GetXaxis()->SetTitle("x [#mum]");
  det.DM->GetYaxis()->SetTitle("y [#mum]");
  det.DM->GetZaxis()->SetTitle("z [#mum]");

  ////////////////////////////////////////////////////////////////////////////////////////    
  //INITIALIZE SPACE CHARGE HISTOGRAM                                                   //
  // NeffH is one of the two methods (the other is a TF3 function) which defines the    //
  // distribution of the space charge in the sensors - unlike a full TCAD simulators    //
  // which calculate the Neff from the properties of traps, the distribution in KDetSim //
  // should be given                                                                    //
  ////////////////////////////////////////////////////////////////////////////////////////  

  det.NeffH=new TH3F("Neff","Neff",det.nx,0,det.nx,det.ny,0,det.ny,det.nz,0,1);
  det.NeffH->GetXaxis()->SetTitle("x [#mum]");
  det.NeffH->GetYaxis()->SetTitle("y [#mum]");
  det.NeffH->GetZaxis()->SetTitle("z [#mum]");
  

  ////////////////////////////////////////////////////////////////////////////////////////    
  //SETTING UP THE ELECTRODES                                                           //
  ////////////////////////////////////////////////////////////////////////////////////////  

  //Back side
  Float_t BackPos[3]={stripPitch*3/2,1,0.5};    //Define the metallized back - the center of the electrode
  Float_t BackSiz[3]={stripPitch*3/2,0.99,0.5}; //Define the metallized back - size of metallizaion 1 micron thick over the back
  det.ElRectangle(BackPos,BackSiz,2,0);         //Make it happen
  
  // Horizonal Strips                                   
  Float_t StripPos[3]={10,substrateThickness-1,0.5}; //Define the strip position - the center of the strips x will change y,z are fixed for all of them
  Float_t StripSiz[3]={stripSize,0.5,0.5};           //Define the 3 strips - the size of the strips; y and z size are 1 um
  Float_t StripX[]={stripPitch/2, 3*stripPitch/2, 5*stripPitch/2};
  Float_t StripS[]={stripSize,stripSize,stripSize};
  
  // Vertical Bars are defined in the same way as before
  Float_t VerticalPos[3]={1,substrateThickness+(BOXthickness + SOIthickness)/2,0.5};
  Float_t VerticalSiz[3]={1,(BOXthickness + SOIthickness)/2,0.5};
  Float_t VerticalX[]={stripPitch/2, 3*stripPitch/2, 5*stripPitch/2};
  Float_t VerticalS[]={1,1,1};


  // The definition of the strips is put into the EG histogram (ElRectangle) 
  for(int i=0; i<3; i++)
  {
    StripPos[0]=StripX[i]; 
    StripSiz[0]=StripS[i]; 
    VerticalPos[0]=VerticalX[i]; 
    VerticalSiz[0]=VerticalS[i]; 
    
    if(i!=10) 
    {
      det.ElRectangle(StripPos,StripSiz,16385,0);  
      // when you define the strip the third argument defines what is the role of the srips; 
      // bit 1= GND, bit2=at high voltage, bit 14 (16384)=electrode for which the Ramo field is calculated 
      // therfore 16385 means that the strip is at ground and it is the strip for which ramo field is calculated 
      // (bit1 & bit14) 
      // The 4th argument is material = leave a default 0
      det.ElRectangle(VerticalPos,VerticalSiz,16385,0);  
    }
    else 
    {
      det.ElRectangle(StripPos,StripSiz,1,0);  
      det.ElRectangle(VerticalPos,VerticalSiz,1,0);  
    }
  }

  // Logic = definition of the logic block = geometry definition
  Float_t LogicPos[4]={80,substrateThickness+1+BOXthickness + SOIthickness/2,0.5};
  Float_t LogicSiz[4]={10,SOIthickness/2-1,0.5};
  Float_t LogicX[4]={stripPitch/4,stripPitch,2*stripPitch,11*stripPitch/4};
  Float_t LogicS[4]={1*stripSize-2,2*stripSize-2,2*stripSize-2,1*stripSize-2};

  for(int i=0; i<4; i++)
  {
    LogicPos[0]=LogicX[i]; 
    LogicSiz[0]=LogicS[i]; 

    //define it in EG histogram
    if (logicConnected) det.ElRectangle(LogicPos,LogicSiz,det->SetElecVolt(0),0);  
    // bit 15 defines the second high voltage - voltage2 variable
    // (that is logic at different votlage than GND)
  }



  ////////////////////////////////////////////////////////////////////////////////////////    
  //SETTING UP THE MATERIAL AND NEFF                                                    //
  ////////////////////////////////////////////////////////////////////////////////////////  
  
  // det->DM->SetBinContent(x,y,z,materialType); materialType: 0: Si, 2: SiO2, 100: Al 
  for(int j=0;j<=det.ny;j++)
    for(int i=0;i<=det.nx;i++) 
    {
      // if outside detector = air
      // space charge outside = 0
      //      	if(i>=200 && j>=100 && j<=300)

      if(j<=substrateThickness)
      {
	det->DM->SetBinContent(i,j,1,0); 
	det->NeffH->SetBinContent(i,j,1,-20);  //(-20 is negative space charge of concentration 2e13)
      }
      
      if((j>substrateThickness && j<substrateThickness+BOXthickness) || j==2 )
      {
	det->DM->SetBinContent(i,j,1,2); 
	det->NeffH->SetBinContent(i,j,1,oxideCharge); 
      }
      
      if(j==1)  det->DM->SetBinContent(i,j,1,100); 

      // Negative space charge under BOX
      if((j <= substrateThickness) && (j >= substrateThickness))
	if(EG->GetBinContent(i,substrateThickness,1) == 0)
	{
	  det->NeffH->SetBinContent(i,j,1,accumulationCharge);
	}
    }

  ////////////////////////////////////////////////////////////////////////////////////////    
  //SET UP DETECTOR VOLTAGES AND CALCULATE FIELD                                        //
  ////////////////////////////////////////////////////////////////////////////////////////  

  det.Voltage = voltage1;
  det.Voltages->Set(1);
  det.Voltages[0] = voltage2;
  det.Mobility=1;   // select mobility model; 1=Canali
  det.diff=1;       // switch on diffusion of free carries

  det->SetBoundaryConditions(); // calculate  boundary condition for the simulated volume
  det->CalField(0);             // calculate  electric field
  det->CalField(1);             // calculate  ramo (weighting) field
    
  det->SetEntryPoint(entryPointX,entryPointY,0); // define entry point of the track 
  det->SetExitPoint(entryPointX,entryPointY+1,1); // define exit point of teh track 
  det->MipIR(100);                                // split the track into 100 segments = charge buckets and perform the simulation of the drift
  
  TCanvas *c1 = new TCanvas("c1","c1",1); 
     
  det->sum.Draw();                                // draw the induced current
  det->neg.Draw("SAME");                          // component from the drift of electrons
  det->pos.Draw("SAME");                          // component from the drift of holes
                                                  // The integral of these histograms gives the charge in 
                                                  // terms of number of charge buckets - so for 100 buckets you should 
                                                  // get in ideal detector 100 (+/- depending on the polarity of the pulses)

  TCanvas *c2 = new TCanvas("c2","c2", 0,600,700,500);
  det.ShowMipIR(100);                             // open a new canvas and show the drift paths
  

}
