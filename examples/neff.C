{
  TF1 * f1 = new TF1( "f1", "[0]+([1]+([2]+([3]+([4]+[5]*x)*x)*x)*x)*x", -1, 1 );
  f1->SetParameter( 0,  0 );
  f1->SetParameter( 1, -0.5 );
  f1->SetParameter( 2,  0.0 );
  f1->SetParameter( 3, -2.0 );
  f1->SetParameter( 4,  0.0 );
  f1->SetParameter( 5,  0.0 );
  f1->Draw("l");
}

