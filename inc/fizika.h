#ifndef SIMUL_FIZIKA
static double TofC   = 1.60217733e-4/3.6e-6;    // From eV to fC
static double e_0    = 1.60217733e-19;           //As
static double SiZ    =  14.00; 
static double SiA    =  28.08;
static double SiRho  =   2.33;                 // Si intrinsic resistivity  
static double Si_a   =   0.1492;
static double Si_m   =   3.25;
static double Si_X1  =   2.87;
static double Si_X0  =   0.2014;
static double Si_C   =  -4.44;
static double ME     =   0.51099906;           // MeV
static double Si_mue =1350.0e-4;                  // cm^2/(Vs)
static double Si_muh = 480.0e-4;                  // cm^2/(Vs)
static double Kboltz = 8.617385e-5;           // eV/K
static double Q_e    = 1.0e-6;
static double perm0  = 8.854187817e-12;        // F/m
static double perm   = 11.7;
static double permDi = 5.7;                  
static double Nd     = 1e20; // atoms/m^3 
static double mh     = ME*0.558;
static double me     = ME*1.08;             //ev/C^2   
static double SIgap  = 1.12;                     //eV
static double Clight = 299792458;              //m/s
static double hbar   = 6.5821220e-22;            //MeV s
static double hbarc  = 197.327053;                //MeV fm
static double Ro=100; //ohm m 

//#define MAXPOINT 3001
#define MAXPOINT 10001
struct sled 
{
float pos[MAXPOINT];
float neg[MAXPOINT];
  float spos;  //sum of positive
  float sneg;  //sum of negative
};
/*struct sled 
{
double pos[MAXPOINT];
double neg[MAXPOINT];
  double spos;  //sum of positive
  double sneg;  //sum of negative
};*/
struct diffusion
{
int index;
float cut;
float diffch;
struct sled xc;
struct sled yc;
int ds;
struct diffusion *prev;
struct diffusion *next;
};

struct segdrift
{
int index;
float totch;
float sp[3];
  //float lp[3];
  //float rp[3];
struct sled xc;
struct sled yc;
struct sled charge;
struct sled time;
int dstrip; //drift strip
struct segdrift *prev;
struct segdrift *next;
struct diffusion *difftrack;
};


#define SIMUL_FIZIKA
#endif








