
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l mob.C

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

  c1.Print( "mob.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.65, 0.4, 0.95, 0.7 );

  TFile * f1 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400_mob1.root" );
  qvsy->SetTitle( "KDetSim edge-on  8 #upoint10^{15 }n/cm^{2}  400V " );
  qvsy->SetStats(0);
  qvsy->SetMaximum(1);
  qvsy->GetXaxis()->SetRangeUser(-75,75);
  qvsy->SetLineColor(1);
  qvsy->SetNdivisions(-506);
  qvsy->Draw();
  lgnd->AddEntry( qvsy, " Canali", "l" );

  TFile * f2 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400_mob2.root" );
  qvsy->SetStats(0);
  qvsy->SetLineColor(8);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " model 2", "l" );

  TFile * f3 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400_mob3.root" );
  qvsy->SetStats(0);
  qvsy->SetLineColor(2);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " Scharf ", "l" );

  TFile * f4 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400_mob4.root" );
  qvsy->SetStats(0);
  qvsy->SetLineColor(4);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " Jacoboni ", "l" );

  lgnd->Draw( "same" );

  c1.Print( "mob.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "mob.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf mob.ps" );
  ierr = system( "rm -f  mob.ps" );
  cout << "evince mob.pdf" << endl;

}
