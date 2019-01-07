
// from KField.cxx

//------------------------------------------------------------------------------
Double_t Mobility( Double_t E, Double_t T, Double_t Charg, Double_t Neff, Int_t which )
{
  Double_t lfm=0,hfm=0;
  Double_t vsatn,vsatp,vsat;
  Double_t betap,betan;
  Double_t alpha;

  switch( which ) // 1 is default: Canali
    {
    case 0:
      alpha=0.72*TMath::Power(T/300,0.065);
      if( Charg > 0 ) {
	Double_t ulp=460*TMath::Power(T/300,-2.18);
	Double_t uminp=45*TMath::Power(T/300,-0.45);
	Double_t Crefp=2.23e17*TMath::Power(T/300,3.2);
	betap=1;
	vsatp=9.05e6*TMath::Sqrt(TMath::TanH(312/T));
	lfm=uminp+(ulp-uminp)/(1+TMath::Power(Neff/Crefp,alpha));
	hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatp,betap),1/betap));
      }
      else {
	Double_t uln=1430*TMath::Power(T/300,-2);
	Double_t uminn=80*TMath::Power(T/300,-0.45);
	Double_t Crefn=1.12e17*TMath::Power(T/300,3.2);
	betan=2;
	vsatn=1.45e7*TMath::Sqrt(TMath::TanH(155/T));
	lfm=uminn+(uln-uminn)/(1+TMath::Power(Neff/Crefn,alpha));
	hfm=2*lfm/(1+TMath::Power(1+TMath::Power(2*lfm*E/vsatn,betan),1/betan));
      }
      break;
    case 1:
      //printf("%e ",par[0]);
      if( Charg > 0 ) {
	lfm=8.54e5*TMath::Power(T,-1.075)*TMath::Exp(1-T/124.);
	vsatp=1.445e7*TMath::Exp(-T/435.9);
	betap=2.49*TMath::Exp(-T/270.3);
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
      }
      else {
	lfm=2.712e8*TMath::Power(T,-2.133);
	vsatn=1.586e7*TMath::Exp(-T/723.6);
	betan=-8.262e-8*TMath::Power(T,3)+6.817e-5*TMath::Power(T,2)-1.847e-2*T+2.429;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
      }
      break;
    case 2:   // WF2
      if( Charg > 0 ) {
	lfm=480;
	vsatp=9.5e6;
	betap=1;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,1/betap),betap);
      }
      else {
	lfm=1350;
	vsatn=1.1e7;
	betan=0.5;
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,1/betan),betan);
      }
      break;
    case 3: // Klanner Scharf
      Double_t bb,cc,E0;
      if( Charg > 0 ) {
	E0=2970*TMath::Power(T/300,5.63);
	bb=9.57e-8*TMath::Power(T/300,-0.155);
	cc=-3.24e-13;
	lfm=457*TMath::Power(T/300,-2.80);
	if(E>E0) hfm=1./(1/lfm+bb*(E-E0)+cc*TMath::Power(E-E0,2)); else hfm=lfm;
      }
      else {
	E0=2970*TMath::Power(T/300,5.63);
	lfm=1430*TMath::Power(T/300,-1.99);
	vsatn=1.05e7*TMath::Power(T/300,-3.02);
	if(E>E0) hfm=1./(1/lfm+1/vsatn*(E-E0)); else hfm=lfm;
      }
      break;
    case 4:   // Jacoboni
      if( Charg > 0 ) {
	lfm = 474 * TMath::Power(T/300., -2.619);
	vsatp = 0.940e7  * TMath::Power(T/300., -0.226);
	betap = 1.181 * TMath::Power(T/300., 0.633 ); // <100> orientation
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatp,betap),1/betap);
      }
      else {
	lfm = 1440*  TMath::Power(T/300., -2.260);
	vsatn = 1.054e7  *  TMath::Power(T/300., -0.602);
	betan = 0.992 *  TMath::Power(T/300., 0.572); // <100> orientation
	hfm=lfm/TMath::Power(1+TMath::Power(lfm*E/vsatn,betan),1/betan);
      }
      break;
    case 9:
      if( Charg > 0 )
	hfm=0;
      else
	hfm=0;
      break;
    case 10: //Diamond parametrization
      if( Charg > 0 ) {
	lfm=2064;
	vsat=14.1e6;
      }
      else {
	lfm=1714;
	vsat=9.6e6;
      };
      hfm = lfm / ( 1+(lfm*E)/vsat );
      break;
    }
  return hfm;

} // Mobility

//------------------------------------------------------------------------------
Double_t DriftVelocity( Double_t E, Double_t Charg, Double_t T, Double_t Neff, Int_t which )
{
  E *= 1e4; // [V/cm]
  return Mobility( E, T, Charg, Neff, which ) * E;
}

//------------------------------------------------------------------------------
void mob()
{
  double E = 10; // [V/um]
  double N = 400; // [1/um^3]

  for( int m = 1; m <= 4; ++m ) { // 1=Canali, 3=Scharf, 4=Jacoboni 
    cout << m;

    for( double T = 253; T < 300; T += 20 ) { // [K]
      double v = DriftVelocity( fabs(E), -1, T, N, m ) * 1e-5; // [um/ns]
      cout << "  " << v;
    }
    cout  << endl;
  }
  /*   253      273      293
    1  109.708  106.279  102.907   Canali
    2  109.637  109.637  109.637
    3  163.232  131.27   107.469   Scharf
    4  107.906  103.404   99.3427  Jacoboni
  */
}
