
// Daniel Pitzl, Dec 2018
// KDetSim edg signal vs depth, scale fluence

// root -l svsf.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  gStyle->SetHistFillStyle(0);

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();
 
  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //               topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "svsf.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.60, 0.7, 0.95, 0.9 );

  TFile * fs8 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V800.root" );
  svsy->SetTitle( "KDetSim edge-on" );
  svsy->GetXaxis()->SetTitle( "smeared height [mm]" );
  svsy->GetYaxis()->SetTitle( "readout signal [relative]" );
  svsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  svsy->SetStats(0);
  svsy->SetLineColor(1);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, "   8 #upoint10^{15 }n/cm^{2} ", "l" );

  TFile * fs16 = TFile::Open( "edg3_px50_F16_lam10_e4_h4_V800.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(4);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, " 16 #upoint10^{15 }n/cm^{2} ", "l" );

  TFile * fs32 = TFile::Open( "edg3_px50_F32_lam10_e4_h4_V800.root" );
  svsy->SetStats(0);
  svsy->SetLineColor(2);
  svsy->Draw("same");
  lgnd->AddEntry( svsy, " 32 #upoint10^{15 }n/cm^{2} ", "l" );

  lgnd->Draw( "same" );

  c1.Print( "svsf.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "svsf.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf svsf.ps" );
  ierr = system( "rm -f  svsf.ps" );
  cout << "evince svsf.pdf" << endl;

}
