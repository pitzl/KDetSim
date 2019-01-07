
// Daniel Pitzl, Dec 2018
// KDetSim edg

// root -l s8-800.C

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

  c1.Print( "s8-800.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.65, 0.70, 0.95, 0.9 );

  TFile * fdata = TFile::Open( "/home/pitzl/eudaq/tele-scope/scopedf34437.root" );
  colpvsy->SetStats(0);
  colpvsy->SetTitle( "8 #upoint10^{15 }n/cm^{2} edge-on 800 V" );
  colpvsy->SetLineColor(1);
  colpvsy->GetXaxis()->SetRangeUser(-0.1,0.1);
  colpvsy->Draw();
  lgnd->AddEntry( colpvsy, " Data 800 V ", "l" );
  int i0 = colpvsy->FindBin(-0.0525);
  double ph0 = colpvsy->GetBinContent(i0);
  cout << "ph0 " << ph0 << endl;

  TH1F * data = (TH1F*)colpvsy->Clone();

  // smeared sim:

  //TFile * fsim = TFile::Open( "edg3_px50_F8_mode1_e3_h3_V800.root" );

  //TFile * fsim = TFile::Open( "edg3_px25_F8_mode3_e4_h1_V800.root" );
  //TFile * fsim = TFile::Open( "edg3_px50_F8_mode3_e4_h2_V800.root" );
  //TFile * fsim = TFile::Open( "edg3_px50_F8_mode3_e4_h3_V800_smr8.root" ); // best
  //TFile * fsim = TFile::Open( "edg3_px40_F8_mode3_e4_h3_V800.root" ); // shallow overshoot
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V800.root" ); // wrong slope
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e4_h3_V800.root" ); // slightly too steep
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e4_h2_V800_smr10.root" ); // good
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e4_h4_V800_smr11.root" ); // under
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e5_h4_V800_smr11.root" ); // over
  //TFile * fsim = TFile::Open( "edg3_px50_F8_lam10_e4p5_h4_V800_smr11.root" ); // best
  TFile * fsim = TFile::Open( "edg3_px50_F8_dlt10_e4_h4_V800.root" ); // mean dy 1.68
  //TFile * fsim = TFile::Open( "edg3_px50_F8_dlt11_e4_h4_V800.root" ); // mean dy 2.84
  //TFile * fsim = TFile::Open( "edg3_px50_F8_dlt10_e6_h6_V800.root" ); // mob4 mean dy 2.00

  //TFile * fsim = TFile::Open( "edg9_px25_F8_mode3_e4_h1_V800.root" ); // 
  //TFile * fsim = TFile::Open( "edg9_px50_F8_mode3_e4_h1_V800.root" ); // rounder, too low
  //TFile * fsim = TFile::Open( "edg9_px50_F8_mode3_e4_h2_V800_smr10.root" ); // 
  //TFile * fsim = TFile::Open( "edg9_px50_F8_lam10_e4_h3_V800.root" ); // overshoot

  // pz 152 shift left:
  //TFile * fsim = TFile::Open( "edg15_px25_F8_dlt10_e3_h2_V800_smr11.root" ); // mean dy 1.30

  //TFile * fsim = TFile::Open( "edg15_px25_F8_dlt11_e2_h2_V800.root" ); // mdy 1.31
  //TFile * fsim = TFile::Open( "edg15_px25_F8_dlt11_e2_h1_V800.root" ); // mdy 1.06

  svsy->SetStats(0);
  svsy->SetLineColor(617); // Magenta
  svsy->SetMarkerColor(617);
  int k0 = svsy->FindBin(-0.0525);
  double qh0 = svsy->GetBinContent(k0);
  cout << "qh0 " << qh0 << endl;
  svsy->Scale(ph0/qh0);
  svsy->Draw( "same" );
  lgnd->AddEntry( svsy, " KDetSim", "l" );

  lgnd->Draw( "same" );

  c1.Print( "s8-800.ps" );

  // data have 5 um y bins

  double sumdy = 0;
  double chisq = 0;
  int npnt = 0;

  for( int ii = 1; ii <= data->GetNbinsX(); ++ii ) {

    double y = data->GetBinCenter(ii);
    if( y < -0.08 ) continue;
    if( y >  0.08 ) break;

    double q = data->GetBinContent(ii);
    double e = data->GetBinError(ii);

    int jj = svsy->FindBin(y);
    double s = svsy->GetBinContent(jj);

    cout << ii
	 << "  " << y*1E3
	 << "  " << svsy->GetBinCenter(jj)*1E3
	 << "  " << q
	 << "  " << s
	 << "  " << e
	 << "  " << s-q
	 << endl;

    sumdy += fabs(s-q);
    chisq += pow( (s-q)/e, 2 );
    ++npnt;

  }

  cout
    << "mean dy " << sumdy/npnt << endl
    << "mean sy " << sqrt(chisq/npnt) << endl
    << "norm ph " << ph0/qh0
    << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "s8-800.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf s8-800.ps" );
  ierr = system( "rm -f  s8-800.ps" );
  cout << "evince s8-800.pdf" << endl;

}
