
// Daniel Pitzl, Dec 2018
// KDetSim vert weighting field

// root -l wvsy.C

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

  c1.Print( "wvsy.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.72, 0.60, 0.95, 0.9 );

  TFile * f1 = TFile::Open( "vert1_lam10_V800.root" );
  wvsy->SetTitle( "KDetSim vertical" );
  wvsy->GetYaxis()->SetTitle( "central pixel weighting potential" );
  wvsy->SetStats(0);
  wvsy->SetMaximum(1);
  wvsy->GetXaxis()->SetRangeUser(-75,75);
  wvsy->SetNdivisions(-506);
  wvsy->SetLineColor(1);
  wvsy->Draw();
  lgnd->AddEntry( wvsy, "   1 pixel", "l" );

  TFile * f3 = TFile::Open( "vert3_lam10_V800.root" );
  wvsy->SetStats(0);
  wvsy->SetLineColor(2);
  wvsy->Draw( "same" );
  lgnd->AddEntry( wvsy, "   3 pixels", "l" );

  TFile * f9 = TFile::Open( "vert9_lam10_V800.root" );
  wvsy->SetStats(0);
  wvsy->SetLineColor(4);
  wvsy->Draw( "same" );
  lgnd->AddEntry( wvsy, "   9 pixels ", "l" );

  TFile * f25 = TFile::Open( "vert25_lam10_V800.root" );
  wvsy->SetStats(0);
  wvsy->SetLineColor(8);
  wvsy->Draw( "same" );
  lgnd->AddEntry( wvsy, " 25 pixels ", "l" );

  lgnd->Draw( "same" );

  c1.Print( "wvsy.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "wvsy.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf wvsy.ps" );
  ierr = system( "rm -f  wvsy.ps" );
  cout << "evince wvsy.pdf" << endl;

}
