#include "KMesh.h"
#include "TMath.h"

ClassImp(KMesh)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KMesh                                                                //
// A mesh generator for the non-equividistant bins                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


Int_t KMesh::GetBins(Int_t Num, Float_t SS, Float_t ES, Float_t *Bins)
{

 Float_t dX=(Max-Min)/Num;
 Float_t dS=(SS-ES)/(Num-1);
 Float_t dSi, SumS=0;
 Int_t Ni=0,k=0; 
 Bins[0]=0;
 N=0;
 printf("Steps=%d, SumS=%f \n",k,SumS); k++;
 for(Int_t i=0;i<Num;i++)
   {
     dSi=SS-i*dS;
     Ni=TMath::Nint(dX/dSi);
     printf("New region : %d %d %f\n",i,Ni,dSi);
     for(Int_t j=0;j<Ni;j++) 
       {
	  SumS+=dSi; Bins[k]=SumS;
	  if(Bins[k]>Max) {Bins[k]=Max; printf("END\n"); break;}
	  printf("Steps=%d, SumS=%f \n",k,Bins[k]);
	 k++;

       }
   }
 if(Bins[k-1]<Max) {Bins[k]=Max;  printf("Steps=%d, SumS=%f \n",k,Bins[k]); };

 return k;

}


Int_t KMesh::GetBins(Int_t size,Float_t *Pos, Float_t *Step, Float_t *Bins)
{
  Float_t N[100];
  Int_t num,NN=0,i,j,k;
  for(i=0;i<size;i++)
    {
      //      if(i>0) N[i]=(Pos[i]-Pos[i-1])/Step[i]; else 
      N[i]=Pos[i]/Step[i];
      if(N[i]!=TMath::Nint(N[i])) {printf("Step size not integer ...\n"); return -1;} 
      else
	NN+=(Int_t)N[i];
    }
  printf("%d \n",NN);
  //  Bins=new Float_t [NN+1];
   
  Bins[0]=0; k=1;
  for(i=0;i<size;i++)
    {
      for(j=0;j<N[i];j++)
	{

	  Bins[k]=Bins[k-1]+Step[i];
	  printf("%d %f\n",k,Bins[k]);
	  k++;
	}
    }

  return NN;
}
