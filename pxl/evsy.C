
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l evsy.C

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

  c1.Print( "evsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.75, 0.4, 0.95, 0.9 );

  // 800 V:

  //TFile * f800 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V800.root" );
  TFile * fs800 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V800.root" ); // mdy 1.06
  evsy->SetTitle( "KDetSim 8E15" );
  evsy->SetStats(0);
  evsy->GetXaxis()->SetRangeUser(-75,75);
  evsy->SetNdivisions(506);
  evsy->SetLineColor(1);
  evsy->Draw();
  lgnd->AddEntry( evsy, " 800 V", "l" );

  // 600 V:

  //TFile * f600 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V600.root" );
  TFile * fs600 = TFile::Open( "edg15_px25_F8_dlt11_e2_h2_V600.root" ); // mdy 1.59
  evsy->SetStats(0);
  evsy->SetLineColor(2);
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 600 V", "l" );

  // 500 V:

  //TFile * f500 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V500.root " );
  TFile * fs500 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V500.root" ); // mean dy 1.09
  evsy->SetStats(0);
  evsy->SetLineColor(4);
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 500 V", "l" );

  // 400 V:

  //TFile * f400 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V400.root " );
  TFile * fs400 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V400.root" ); // mdy 1.41
  evsy->SetStats(0);
  evsy->SetLineColor(8);
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 400 V", "l" );

  // 300 V:

  //TFile * f300 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V300.root" );
  TFile * fs300 = TFile::Open( "edg15_px25_F8_dlt11_e2_h3_V300.root" ); // mdy 1.18
  evsy->SetStats(0);
  evsy->SetLineColor(6); // magenta
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 300 V", "l" );

  // 250 V:

  //TFile * f250 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V250.root" );
  TFile * fs250 = TFile::Open( "edg15_px25_F8_dlt11_e3_h3_V250.root" ); // mean dy 1.52
  evsy->SetStats(0);
  evsy->SetLineColor(7); // cyan
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 250 V", "l" );

  // 200 V:

  //TFile * f200 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V200.root" );
  TFile * fs200 = TFile::Open( "edg15_px25_F8_dlt11_e3_h4_V200.root" ); // mdy 1.02
  evsy->SetStats(0);
  evsy->SetLineColor(94); // orange
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 200 V", "l" );

  // 150 V:

  //TFile * f150 = TFile::Open( "edg3_px25_F8_mode3_e4_h2_V150.root" );
  TFile * fs150 = TFile::Open( "edg15_px25_F8_dlt11_e4_h3_V150.root" ); // mdy 1.10
  evsy->SetStats(0);
  evsy->SetLineColor(11); // grey
  evsy->Draw( "same" );
  lgnd->AddEntry( evsy, " 150 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "evsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "evsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf evsy.ps" );
  ierr = system( "rm -f  evsy.ps" );
  cout << "evince evsy.pdf" << endl;

}
