
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l rvsy.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

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

  c1.Print( "rvsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.82, 0.5, 0.95, 0.9 );

  // 800 V:

  TFile * f800 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V800.root" );
  rvsy->SetTitle( "KDetSim edge-on" );
  rvsy->SetStats(0);
  rvsy->GetXaxis()->SetRangeUser(-75,75);
  rvsy->SetNdivisions(-506);
  rvsy->SetLineColor(1);
  rvsy->Draw();
  lgnd->AddEntry( rvsy, " 800 V", "l" );

  // 600 V:

  TFile * f600 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V600.root" );
  rvsy->SetStats(0);
  rvsy->SetLineColor(2);
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 600 V", "l" );

  // 500 V:

  TFile * f500 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V500.root " );
  rvsy->SetStats(0);
  rvsy->SetLineColor(4);
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 500 V", "l" );

  // 400 V:

  TFile * f400 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V400.root " );
  rvsy->SetStats(0);
  rvsy->SetLineColor(8);
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 400 V", "l" );

  // 300 V:

  TFile * f300 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V300.root" );
  rvsy->SetStats(0);
  rvsy->SetLineColor(6); // magenta
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 300 V", "l" );

  // 250 V:

  TFile * f250 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V250.root" );
  rvsy->SetStats(0);
  rvsy->SetLineColor(7); // cyan
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 250 V", "l" );

  // 200 V:

  TFile * f200 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V200.root" );
  rvsy->SetStats(0);
  rvsy->SetLineColor(94); // orange
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 200 V", "l" );

  // 150 V:

  TFile * f150 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V150.root" );
  rvsy->SetStats(0);
  rvsy->SetLineColor(9); // sky blue
  rvsy->SetLineColor(28); // brown
  rvsy->Draw( "same" );
  lgnd->AddEntry( rvsy, " 150 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "rvsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "rvsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf rvsy.ps" );
  ierr = system( "rm -f  rvsy.ps" );
  cout << "evince rvsy.pdf" << endl;

}
