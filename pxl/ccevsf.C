
// Daniel Pitzl, Dec 2018
// KDetSim vertical charge collecion efficiency

// root -l ccevsf.C

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

  c1.Print( "ccevsf.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.65, 0.7, 0.95, 0.9 );

  TFile * fs12 = TFile::Open( "vert9_lam10_V1200.root" );
  qvsf->SetTitle( "exponential KDetSim vertical" );
  qvsf->SetStats(0);
  qvsf->SetMarkerStyle(20);
  qvsf->SetMarkerColor(2);
  qvsf->Draw( "psame" );
  lgnd->AddEntry( qvsf, " 1200 V", "p" );

  TFile * fs8 = TFile::Open( "vert9_lam10_V800.root" );
  qvsf->SetStats(0);
  qvsf->SetMarkerStyle(20);
  qvsf->SetMarkerColor(1);
  qvsf->Draw( "psame" );
  lgnd->AddEntry( qvsf, "   800 V", "p" );

  TFile * fs4 = TFile::Open( "vert9_lam10_V400.root" );
  qvsf->SetStats(0);
  qvsf->SetMarkerStyle(20);
  qvsf->SetMarkerColor(4);
  qvsf->Draw("psame");
  lgnd->AddEntry( qvsf, "   400 V", "p" );

  lgnd->Draw( "same" );

  c1.Print( "ccevsf.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "ccevsf.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf ccevsf.ps" );
  ierr = system( "rm -f  ccevsf.ps" );
  cout << "evince ccevsf.pdf" << endl;

}
