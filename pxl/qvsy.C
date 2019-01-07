
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l qvsy.C

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

  c1.Print( "qvsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.80, 0.5, 0.95, 0.9 );

  // 800 V:

  //TFile * fs800 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V800.root" ); // best
  TFile * fs800 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V800.root" ); // mdy 1.06
  qvsy->SetStats(0);
  qvsy->SetTitle( "KDetSim 8E15" );
  qvsy->SetLineColor(1);
  qvsy->GetXaxis()->SetRangeUser(-75,75);
  qvsy->SetNdivisions(-506);
  qvsy->Draw();
  lgnd->AddEntry( qvsy, " 800 V", "l" );

  // 600 V:

  //TFile * fs600 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V600.root" );
  TFile * fs600 = TFile::Open( "edg15_px25_F8_dlt11_e2_h2_V600.root" ); // mdy 1.59
  qvsy->SetStats(0);
  qvsy->SetLineColor(2);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 600 V", "l" );

  // 500 V:

  //TFile * fs500 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V500.root" );
  TFile * fs500 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V500.root" ); // mean dy 1.09
  qvsy->SetStats(0);
  qvsy->SetLineColor(4);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 500 V", "l" );

  // 400 V:

  //TFile * fs400 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V400.root" );
  TFile * fs400 = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V400.root" ); // mdy 1.41
  qvsy->SetStats(0);
  qvsy->SetLineColor(8);
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 400 V", "l" );

  // 300 V:

  //TFile * fs300 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V300.root" );
  TFile * fs300 = TFile::Open( "edg15_px25_F8_dlt11_e2_h3_V300.root" ); // mdy 1.18
  qvsy->SetStats(0);
  qvsy->SetLineColor(6); // magenta
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 300 V", "l" );

  // 250 V:

  //TFile * fs250 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V250.root" );
  TFile * fs250 = TFile::Open( "edg15_px25_F8_dlt11_e3_h3_V250.root" ); // mean dy 1.52
  qvsy->SetStats(0);
  qvsy->SetLineColor(7); // cyan
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 250 V", "l" );

  // 200 V:

  //TFile * fs200 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V200.root" );
  TFile * fs200 = TFile::Open( "edg15_px25_F8_dlt11_e3_h4_V200.root" ); // mdy 1.02
  qvsy->SetStats(0);
  qvsy->SetLineColor(94); // orange
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 200 V", "l" );

  // 150 V:

  //TFile * fs150 = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V150.root" );
  TFile * fs150 = TFile::Open( "edg15_px25_F8_dlt11_e4_h3_V150.root" ); // mdy 1.10
  qvsy->SetStats(0);
  qvsy->SetLineColor(11); // grey
  qvsy->Draw( "same" );
  lgnd->AddEntry( qvsy, " 150 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "qvsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "qvsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf qvsy.ps" );
  ierr = system( "rm -f  qvsy.ps" );
  cout << "evince qvsy.pdf" << endl;

}
