
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l taue.C

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

  c1.Print( "taue.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.72, 0.60, 0.95, 0.9 );

  TFile * f200 = TFile::Open( "edg3_px25_F8_mode3_e8_h2_V800.root" );
  qvsy->SetTitle( "KDetSim edge-on" );
  qvsy->SetStats(0);
  qvsy->SetMaximum(1);
  qvsy->GetXaxis()->SetRangeUser(-75,75);
  qvsy->SetNdivisions(-506);
  qvsy->SetLineColor(1);
  qvsy->Draw();
  lgnd->AddEntry( qvsy, " #tau_{e}  1 ns", "l" );

  TFile * f20 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V800.root" );
  qvsy->SetStats(0);
  qvsy->SetLineColor(2);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " #tau_{e}  0.5 ns", "l" );

  TFile * f2 = TFile::Open( "edg3_px25_F8_mode3_e2_h2_V800.root " );
  qvsy->SetStats(0);
  qvsy->SetLineColor(4);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " #tau_{e}  0.25 ns ", "l" );

  lgnd->Draw( "same" );

  c1.Print( "taue.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "taue.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf taue.ps" );
  ierr = system( "rm -f  taue.ps" );
  cout << "evince taue.pdf" << endl;

}
