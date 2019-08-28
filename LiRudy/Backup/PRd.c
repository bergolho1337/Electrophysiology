////////////////////////////////////////////////////////////////////////
//     The Pan-Rudy Dynamic (PRd) Model of the Canine Purkinje Cell
//                           as published in 
// "A model of Canine Purkinje Electrophysiology and Ca Cycling: Rate
// Dependence, Triggered Activity and Comparison to Ventricular Myocyte"
//         Li P and Rudy Y. Circulation Research. 2011;109;71-79
//      (http://circres.ahajournals.org/cgi/content/full/109/1/71)
//             Email: panli@seas.wustl.edu   rudy@wustl.edu
//            Writted by Pan Li, Last Modified: 24 June 2011
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define	 beats 10       
#define	 BCL   1000	  
#define	 S2    1000

// GEOMETRY
double ageo,acap,vcell,vmyo,vnsr,vjsr,vsss,vcsr,vssl,vmito;
const double pi = 3.14;			
const double radius = 0.00175;	
const double length = 0.0164;	
const double rcg = 1.54;	

// PHYSICAL CONSTANTS
const double frdy = 96485;		      
const double R = 8314;			  
const double temp = 310;		    

const double nao = 140;			
const double cao = 1.8;			
const double ko  = 5.4;			
const double clo = 100;		

const double zna = 1;			
const double zk  = 1;			
const double zcl = -1;		
const double zca = 2;			
const double ganai = 0.75;		
const double ganao = 0.75;		
const double gaki  = 0.75;	
const double gako  = 0.75;	
const double gacai = 1.0;		
const double gacao = 0.341;	

// REVERSAL POTENTIALS
double ena,ek,eca;

// TIME
double t,dt,tmax,T;				
const double dtmin = 0.001;		     
const double dtmed = 0.005;		
const double dtmax = 0.1;		

// VOLTAGE
double v,dvdt,dvdtclock;			
const double dvdtthresh = 1;		

// STIMULUS CURRENT
const double stimdur = 0.5;	        
double istim = -80;			        
double tstim,stimtime;			   
int stimcount;				        

// TOTAL TRANSMEMBRANE CURRENTS
double icatot,iktot,inatot,icltot,itot;					

// MEMBRANE IONIC CURRENTS
double ina;				
const double gna = 18;	
double m,ma,mb,mtau,mss,h,ha,hb,htau,hss,j,ja,jb,jtau,jss;			

double inal;			
const double gnal2 = 0.052;	
const double gnal3 = 0.018;	
double ml,alphaml,betaml,mltau,mlss,hl,hltau,hlss;	
double jltau, jlss, jl, jl3tau, jl3ss, jl3;
double inal2, inal3;
double i3tau, i3ss, i3,Rtau, Rss, Ri,ml3tau, ml3ss, hl3tau, hl3, ml3, hl3ss,ireltau, REL;				

double ical,ibarca;				
const double pca = 1.9926e-4;	
double d,dss,dtau,dpower,powss;				
const double powtau = 10;	
double f,fss,ftau,f2,f2ss,f2tau,fca,fcass,fcatau,fca2,fca2tau,fca2ss;			

double icat;    		         
const double gcat = 0.07875;    	
double taub,taug,bss,gss,b,g;       	

double icab;				

double ito1;			
double itos, itof;	
const double gtos = 0.1414;	   
const double gtof = 0.042;    
double rto1,a,alphaa,betaa,atau,ass,i,alphai,betai,itau,iss,i2,alphai2,betai2,i2tau,i2ss;				

double ikr,gkr,xr,xrss,xrtau,rkr;				

double iks,gks,eks;				
const double prnak = 0.014;		
double xs1,xs1tau,xs2,xs2tau,xsss;

double ik1,gk1,k1ss;	
	
double yss, ytau, y, ifna, ifk, iftotal;

double inab;				    
const double gnab = 0.0025;     

const double pcab = 3.99e-8;	  
const double pnab = 0.64e-8;	

double inaca;				    
const double inacamax = 2.52;
const double kmcaact = 0.000125;
const double kmnai1 = 12.3;		
const double kmnao = 87.5;		
const double kmcai = 0.0036;	
const double kmcao = 1.3;		
const double nu = 0.35;			
const double ksat = 0.27;		
double allo,num,denommult,denomterm1,denomterm2,deltaE,inacass;				

double inak;	
const double ibarnak = 1.1004;		

double ipca;				
const double ipcabar = 0.0115;		   
const double kmpca = 0.0005;

// CALCIUM FLUXES RATE CONSTANTS
const double tautr1 = 120;
const double tautr2 = 120;	
const double gaptau = 12;
const double sstau = 0.2;	

// CALCIUM FLUXES AND CONCENTRATIONS
double qrel1,qrel2;			
double irelss;	
double IP3 = 0.0001;  
const double k1 = 150000;
const double k1a = 16.5;
const double k0 = 96000;
const double k0a = 9.6;
const double k2 = 1800;
const double k2a = 0.21;
const double tauip3r = 3.7;
double du,u,POip3,qip3;

const double  dqupcamkbar = 0.75;	
double  dqupcamk;			        
const double  dkmplbbar = 0.00017;	
double  dkmplb,qup1,qup2;	
double kmup   = 0.00028;	
const double nsrbar = 15.0;	 

double  qtr1,qtr2;			

double bsss,csqn1,bjsr,cjsr,csqn,bcsr,ccsr;
double casss,dcasss,cassstot,bsr,bsl,b1,c1,d1;				
const double bsrbar = 0.019975;	
const double kmbsr = 0.00087;		
const double bslbar = 0.4777;	
const double kmbsl = 0.0087;	
double cassl,dcassl;	
double qdiff;				               

double cajsr,dcajsr,cajsrtot;	     
const double csqnbar = 2.88;		    
const double kmcsqn = 0.8;		   
double cacsr,dcacsr,cacsrtot;			         
double cansr,dcansr;

double cai,dcai,catotal,cmdn;				     
const double cmdnbar = 0.1125;	
const double kmcmdn = 2.38e-3;	
double trpn;				        
const double trpnbar = 3.15e-2;
const double kmtrpn = 0.5e-3;	
double bmyo,cmyo,dmyo;				   
double trpnbar1 = 3.5e-3;
double cmdnbar1 = 1.25e-2;
double csqnbar1 = 1.2; 			

double caavg;	         

// SODIUM/POTASSIUM FLUXES AND CONCENTRATIONS
double nai,dnai,nasss,dnasss,ki,dki,ksss,dksss;	
double qgap, qdiffna, qgapna, nassl, dnassl;

// CAMKII DYNAMICS
double camkbound,camkactive,camktrap;	
const double camk0 = 0.05;		
const double alphacamk = 0.05;		
const double betacamk = 0.00068;	
const double kmcam = 0.0015;		
const double kmcamk = 0.15;		   
double fca_dtaucamk;			
const double fca_dtaucamkbar = 10.0;

// OUTPUT FILE
FILE *ap;		
int count;			

// FUNCTIONS
void timestep ();			
void comp_revs ();			
void comp_ina ();		
void comp_ical ();		
void comp_icat ();
void comp_ik1 ();			
void comp_iks ();		
void comp_ikp ();			
void comp_icab ();			
void comp_ikb ();
void comp_inab ();
void comp_ikr ();		
void comp_inaca ();			
void comp_inak ();			
void comp_inal ();			
void comp_ipca ();			
void comp_ito1 ();			
void comp_istim ();			
void comp_itot ();			
void comp_qup1 ();			
void comp_qup2 ();	
void comp_qrel1 ();			
void comp_qrel2 ();
void comp_ip3 ();
void comp_qtr1 ();			
void comp_qtr2 ();
void comp_conc ();			
void comp_if ();
void printtofile ();			

int main ()
{
	// OPEN OUTPUT FILE
	ap =  fopen("PRd.txt","w");
	
	// CELL GEOMETRY
	vcell	= 1000*pi*radius*radius*length;
	ageo	= 2*pi*radius*radius + 2*pi*radius*length;
	acap	= rcg*ageo;
	vmyo	= vcell * 0.60;
	vnsr	= vcell * 0.04;
	vmito   = vcell * 0.18;
	vjsr	= vcell * 0.002;
	vcsr	= vcell * 0.008;
	vsss	= vcell * 0.02;
	vssl    = vcell * 0.15;
	
	// INITIAL CONDITIONS
	/*
    v		= -85;
    m		= 0.0;
    h		= 0.9;
    j		= 0.9;
    d		= 0.0;
    f		= 0.9;
    f2		= 0.9;
    fca		= 0.9;
    fca2	= 0.9;
    xs1		= 0.0;
    xs2		= 0.0;
    xr		= 0.0;
    a		= 0.0;
    i		= 0.9;
    i2		= 0.9;
    ml		= 0.0;
    ml3		= 0.0;
    hl		= 0.9;
    hl3		= 0.9;
    jl		= 0.9;
    jl3     = 0.9;
    casss	= 0.0001;
    cajsr	= 1.0;
    cacsr	= 1.0;
    cansr	= 1.0;
    cassl	= 0.0001; 
    nai		= 8.0;
    nassl	= 8.0;
    nasss	= 8.0;
    ki		= 140.0;
    cai		= 0.0001;
    b	    = 0.0;
    g	    = 0.9;
    u       = 0.0;
    y       = 0.0;
    camktrap= 0.0;
    */
    
    // STATE STATE CONDITIONS DURING PACING CL=1000ms
    v		= -84.058830;
    m		= 0.000821;
    h		= 0.995741;
    j		= 0.999872;
    d		= 0.000016;
    f		= 0.999193;
    f2		= 0.988692;
    fca		= 0.965405;
    fca2	= 0.739378;
    xs1		= 0.001114;
    xs2		= 0.042234;
    xr		= 0.069808;
    a		= 0.000119;
    i		= 0.992541;
    i2		= 0.745628;
    ml		= 0.000329;
    ml3		= 0.046538;
    hl		= 0.984170;
    hl3		= 0.853893;
    jl		= 0.912569;
    jl3		= 0.827885;
    casss	= 0.000135;
    cajsr	= 1.510741;
    cacsr	= 1.537577;
    cansr	= 1.538668;
    cassl	= 0.000130; 
    nai		= 11.501546;
    nassl	= 11.501230;
    nasss	= 11.501240;
    ki		= 136.422946;
    cai		= 0.000053;
    b	    = 0.000437;
    g	    = 0.990384;
    u       = 0.535627;
    y       = 0.182859;
    camktrap= 0.010600;
    

    // TIME SETTINGS
	dvdtclock	= 1000;
	tmax		= BCL*(beats)+S2+500;
	t		= 0;
	T       = 0;
	tstim		= 0;
	stimtime	= 1000;
	stimcount 	= -1;
	
	// MAIN LOOP
	while (t<=tmax) 	
	{

		timestep ();
		
		comp_revs ();
		comp_ina ();
		comp_inal ();
		comp_inab ();
        comp_ical ();
		comp_icat ();
		comp_icab ();
		comp_ito1 ();
		comp_ikr ();
		comp_iks ();
		comp_ik1 ();
		comp_inaca ();
		comp_inak ();
		comp_ipca ();
		comp_if ();
		comp_istim ();
		comp_itot ();
		
		comp_ip3 ();
		comp_qrel1 ();
		comp_qrel2 ();
		comp_qup1 ();
		comp_qup2 ();
		comp_qtr1 ();
		comp_qtr2 ();
		
		comp_conc ();
		
		dvdt	= -itot;
		v	   += dvdt*dt;
       	
		printtofile ();
		
		t	+= dt;

	}

	fclose(ap);
	return(1);
} 

void timestep ()
{
	if ((dvdt>dvdtthresh) || (t>(tstim-2)) || (stimtime<(stimdur+2)) || (dvdtclock<5))
		dt = dtmin; 			
	else if (dvdtclock>=5 && dvdtclock<20)
		dt = dtmed; 			
	else
		dt = dtmax; 			
	
	if (dvdt>dvdtthresh) dvdtclock = 0;	
	dvdtclock += dt;			
}

void comp_revs ()
{
    eca	= (R*temp/(zca*frdy))*log(cao/cassl);
	ena	= (R*temp/frdy)*log(nao/nassl);
	ek	= (R*temp/frdy)*log(ko/ki);
}

void comp_ina ()
{ 
    ma	= 0.64*(v+37.13)/(1-exp(-0.1*(v+37.13)));
	mb	= 0.16*exp(-v/11);
	if (v<-40)
	{
		ha = 0.135*exp((70+v)/-6.8);
		hb = 3.56*exp(0.079*v)+310000*exp(0.35*v);
		ja = (-127140*exp(0.2444*v)-0.003474*exp(-0.04391*v))*(v+37.78)/(1+exp(0.311*(v+79.23)));
		jb = 0.1212*exp(-0.01052*v)/(1+exp(-0.1378*(v+40.14)));
	}
	else
	{
		ha = 0.0;
		hb = 1/(0.13*(1+exp((v+10.66)/-11.1)));
		ja = 0.0;
		jb = 0.3*exp(-0.0000002535*v)/(1+exp(-0.1*(v+32)));
	}
	mtau	= 1/(ma+mb);
	htau	= 1/(ha+hb);
	jtau	= 1/(ja+jb);
	mss	= ma*mtau;
	hss	= ha*htau;
	jss	= 1*ja*jtau;
	m	= mss-(mss-m)*exp(-dt/mtau);
	h	= hss-(hss-h)*exp(-dt/htau);
	j	= jss-(jss-j)*exp(-dt/jtau);
	ina	= gna*pow(m,3)*h*j*(v-ena);
}

void comp_inal ()
{
	mltau	= 1/(0.64*(v+37.13)/(1-exp(-0.1*(v+37.13))) + 0.16*exp(-v/11));
	ml3tau  = mltau;
	mlss	= 1/(1+exp(-(v+28)/7));
	ml3ss   = 1/(1+exp(-(v+63)/7));
	hltau   = 162+132/(1+exp(-(v+28)/5.5));
	hl3tau  = 0.5*hltau;
	hlss	= 1/(1+exp((v+28)/12));
	hl3ss	= 1/(1+exp((v+63)/12));
	jltau   = 411;
	jl3tau  = 0.5*jltau;
	jlss	= hlss;
	jl3ss	= hl3ss;
	ml	    = mlss-(mlss-ml)*exp(-dt/mltau);
	ml3     = ml3ss-(ml3ss-ml3)*exp(-dt/ml3tau);
	hl	    = hlss-(hlss-hl)*exp(-dt/hltau);
	hl3     = hl3ss-(hl3ss-hl3)*exp(-dt/hl3tau);
	jl	    = jlss-(jlss-jl)*exp(-dt/jltau);
	jl3     = jl3ss-(jl3ss-jl3)*exp(-dt/jl3tau);
	inal2   = gnal2*ml*hl*jl*(v-ena);
	inal3   = gnal3*ml3*hl3*jl3*(v-ena);
	inal    = inal2 + inal3; 
}

void comp_ical ()
{
	ibarca		= pca*zca*zca*(((v-15)*frdy*frdy)/(R*temp))*((gacai*casss*exp((zca*(v-15)*frdy)/(R*temp))-gacao*cao)/(exp((zca*(v-15)*frdy)/(R*temp))-1));
	dss		    = (1/(1.0+exp(-(v-2.0)/7.8)));
	dtau		= (0.59+0.8*exp(0.052*(v+13))/(1+exp(0.132*(v+13))));
	fss	        = 1/(1.0 + exp((v+16.5)/9.5));
	ftau        = 0.92/(0.125*exp(-(0.058*(v-2.5))*(0.045*(v-2.5)))+0.1);
	f2ss        = fss;
	f2tau       = 0.90/(0.02*exp(-(0.04*(v-18.6))*(0.045*(v-18.6)))+0.005);
	fcass		= 0.3/(1 - ical/0.05) + 0.55/(1.0+casss/0.003)+0.15;
	fcatau		= 10*camkactive/(camkactive+kmcam) + 0.5+1/(1.0+casss/0.003);
	fca2ss		= 1.0/(1.0-ical/0.01);
	fca2tau		= 1*(300.0/(1.0+exp((-ical-0.175)/0.04))+125.0);
	d		    = dss-(dss-d)*exp(-dt/dtau);
	f		    = fss-(fss-f)*exp(-dt/ftau);
	f2		    = f2ss-(f2ss-f2)*exp(-dt/f2tau);
	fca		    = fcass-(fcass-fca)*exp(-dt/fcatau);
	fca2		= fca2ss-(fca2ss-fca2)*exp(-dt/fca2tau);
	ical		= d*f*f2*fca*fca2*ibarca;	
}

void comp_icat ()
{
	bss	    = 1/(1+ exp (-(v+30)/7));
	gss	    = 1/(1+exp((v+61)/5));
	taub	= 1/(1.068*exp((v+16.3)/30)+1.068*exp(-(v+16.3)/30));
	taug    = 1/(0.015*exp(-(v+71.7)/83.3)+0.015*exp((v+71.7)/15.4));
	b	    = bss-(bss-b)*exp(-dt/taub);
	g	    = gss-(gss-g)*exp(-dt/taug);
	icat	= gcat*b*g*(v-eca);
}

void comp_ito1 ()
{
	atau	= 1/(25*exp((v-82)/18)/(1+exp((v-82)/18))+25*exp(-(v+52)/18)/(1+exp(-(v+52)/18)));
	itau	= 2.86+ 1/(exp(-(v+125)/15)*0.1 + 0.1*exp((v+2)/26.5));
	i2tau	= 21.5+ 1/(exp(-(v+138.2)/52)*0.005 + 0.003*exp((v+18)/12.5));
	ass	    = 1/(1+exp(-(v-8.9)/10.3));
	iss	    = 1/(1+exp((v+30)/11));
	i2ss	= iss;
	a	    = ass-(ass-a)*exp(-dt/atau);
	i	    = iss-(iss-i)*exp(-dt/itau);
	i2	    = i2ss-(i2ss-i2)*exp(-dt/i2tau);
	itos    = gtos*a*i*i2*(v-ek);
	itof    = gtof*(v-ek)/(1+exp(-(v-3)/19.8));
	ito1	= itos + itof;
}

void comp_ikr ()
{
	gkr	    = 0.0326*sqrt(ko/5.4);
	xrss	= 1/(1+exp(-(v)/15));
	xrtau   = 400.0/(1.0+exp(v/10.0)) + 100.0;
	rkr	    = 1/(1+exp((v)/35));
	xr	    = xrss-(xrss-xr)*exp(-dt/xrtau);
	ikr	    = gkr*xr*rkr*(v-ek); 
}

void comp_iks ()
{
	eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(ki+prnak*nassl));
	gks	    = 0.053*(1+0.6/(1+pow((0.000038/cassl),1.4)));
	xsss	= 1/(1+exp(-(v-9)/13.7));
	xs1tau	= 200/(exp(-(v+10)/6) + exp((v-62)/55));
	xs2tau	= 1500+ 350/(exp(-(v+10)/4) + exp((v-90)/58));
	xs1	    = xsss-(xsss-xs1)*exp(-dt/xs1tau);
	xs2	    = xsss-(xsss-xs2)*exp(-dt/xs2tau);
	iks	    = gks*xs1*xs2*(v-eks);
}

void comp_ik1 ()
{	
    k1ss      = 1/(1+exp((v+103-(2.9+ko*2.175))/10.15));
	gk1	      = 0.12*sqrt(ko);
	ik1	      = gk1*k1ss*(v-ek);
}

void comp_if ()
{
	yss       = 1/(1+exp((v+87)/9.5));
	ytau      = 2000/(exp(-(v+132)/10) + exp((v+57)/60));
	y         = yss - (yss-y)*exp(-dt/ytau);
	ifna	  = 0.012*y*y*(v-ena);
	ifk       = 0.024*y*y*(v-ek);
	iftotal   = ifna + ifk;
}

void comp_inaca ()
{
	allo		= 1/(1+pow((kmcaact/(1.5*casss)),2));
	num		    = inacamax*(pow(nasss,3)*cao*exp(nu*v*frdy/(R*temp))-pow(nao,3)*1.5*casss*exp((nu-1)*v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v*frdy/(R*temp));
	denomterm1	= kmcao*pow(nasss,3)+pow(kmnao,3)*1.5*casss+pow(kmnai1,3)*cao*(1+1.5*casss/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nasss/kmnai1,3))+pow(nasss,3)*cao+pow(nao,3)*1.5*casss;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inacass		= 0.2*allo*deltaE;
	
	allo		= 1/(1+pow((kmcaact/(1.5*cassl)),2));
	num		    = inacamax*(pow(nassl,3)*cao*exp(nu*v*frdy/(R*temp))-pow(nao,3)*1.5*cassl*exp((nu-1)*v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v*frdy/(R*temp));
	denomterm1	= kmcao*pow(nassl,3)+pow(kmnao,3)*1.5*cassl+pow(kmnai1,3)*cao*(1+1.5*cassl/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nassl/kmnai1,3))+pow(nassl,3)*cao+pow(nao,3)*1.5*cassl;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inaca		= 0.8*allo*deltaE;
}

void comp_inak ()
{
    inak	= ibarnak*(1/(1+exp(-1*(v+92)*frdy/(R*temp))))*pow((nassl/(nassl+2.6)),3)*(ko/(ko+0.8));
}

void comp_ipca ()
{
	ipca	= ipcabar/((kmpca/cassl)+1);
}

void comp_icab ()
{
	icab	= pcab*zca*zca*((v*frdy*frdy)/(R*temp))*((gacai*cassl*exp((zca*v*frdy)/(R*temp))-gacao*cao)/(exp((zca*v*frdy)/(R*temp))-1));
}

void comp_inab ()
{
    inab    = pnab*frdy*((frdy*v)/(R*temp))*(nassl*exp((frdy*v)/(R*temp)) - nao)/(exp((frdy*v)/(R*temp))-1);     
}

void comp_istim () 
{
	stimtime += dt;
	if (t>=tstim)
	{
		stimtime = 0.0;					
		stimcount += 1;					
		if (stimcount < beats-1)  tstim += BCL;		
		else if (stimcount == beats-1) tstim += S2;	
		else tstim = tmax+1;				
		if (stimcount < beats) printf ("S1 Beat %d at time = %.2f ms !\n", stimcount+1, t);
		else if (stimcount == beats) printf ("S2 Beat at time = %.2f ms !\n", t);
	}
}

void comp_itot ()
{
	if (stimtime>=0.0 && stimtime<stimdur)
	{
		icatot	= ical+icat+ipca+icab-2*inaca-2*inacass;
		iktot	= ikr+iks+ik1-2*inak+ito1+ifk+1*istim;
		inatot	= 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
		itot	= icatot+iktot+inatot;
	}
	else
	{
		icatot	= ical+icat+ipca+icab-2*inaca-2*inacass;
		iktot	= ikr+iks+ik1-2*inak+ito1+ifk;
		inatot	= 3*inak+ina+3*inaca+3*inacass+inal+ifna+inab;
		itot	= icatot+iktot+inatot;
	}
}

void comp_ip3 ()
{
    du= dt*(casss*k2*(1-u) - k2a*u);
    u += du;
    POip3 = tauip3r*IP3*casss*(1-u)/((1+IP3*k0/k0a)*(1+casss*k1/k1a));
    qip3 = 10.920*(cajsr-casss)*(POip3);
}

void comp_qrel1 ()
{
	qdiff  = (casss-cassl)/sstau;  
    REL  = -((ical)*acap/(vsss*2.0*frdy) - (qrel1 + qip3)*vjsr/vsss + qdiff);     
    ireltau = 2*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cajsr));
    if (REL > 0)
    {irelss  = 15*(1+1*(1/(1+pow((0.28/camkactive),8))))*REL/(1 + pow((1.0/cajsr),8));}
    else {irelss = 0;}
    qrel1 += dt*((irelss-qrel1)/ireltau);
}

void comp_qrel2 ()
{
	qgap  = (cassl-cai)/gaptau;  
    REL  = (-qup2*vnsr/vmyo + qgap*vssl/vmyo+ (qrel2)*vcsr/vmyo);    
    ireltau = 6*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cacsr));
    if (REL > 0)
    {irelss  = 91*(1+1*(1/(1+pow((0.28/camkactive),8))))*(REL)/(1 + pow((1/cacsr),8));}
    else {irelss = 0;}
    qrel2 += dt*((irelss-qrel2)/ireltau);
}

void comp_qup1 ()
{
    dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	qup1		= 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cassl,1))-0.00105*cansr/nsrbar;
}

void comp_qup2 ()
{
    dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	qup2		= 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cai,1))-0.0042*cansr/nsrbar;
}

void comp_qtr1 ()
{
	qtr1		= (cansr-cajsr)/tautr1;
}

void comp_qtr2 ()
{
	qtr2		= (cansr-cacsr)/tautr2;
}

void comp_conc ()
{
	qdiff       = (casss-cassl)/sstau;  
	qgap        = (cassl-cai)/gaptau;  
    qdiffna     = (nasss-nassl)/sstau;
    qgapna      = (nassl-nai)/gaptau;
    
	dcasss		= dt*(-(ical-2*inacass)*acap/(vsss*2.0*frdy)+(qrel1+qip3)*vjsr/vsss-qdiff);
	bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+casss,2))+(bslbar*kmbsl/pow(kmbsl+casss,2)));
	casss      += bsss*dcasss;
	
	dcassl		= dt*(-(qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(icat+ipca+icab-2*inaca)*acap/(vssl*2.0*frdy));
	trpn        = trpnbar1*(cassl/(cassl+kmtrpn));
	cmdn		= cmdnbar1*(cassl/(cassl+kmcmdn));
	catotal		= trpn+cmdn+dcassl+cassl;
	bmyo		= cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cassl		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;   
 	
	dcajsr		= dt*(qtr1-qrel1-qip3);
	csqn1       = csqnbar1*(cajsr/(cajsr+kmcsqn));
	bjsr        = csqnbar1 - csqn1-cajsr-dcajsr+kmcsqn;
	cjsr        = kmcsqn*(csqn1+cajsr+dcajsr);
	cajsr       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	
	dcacsr		= dt*(qtr2-qrel2);
	csqn        = csqnbar*(cacsr/(cacsr+kmcsqn));
	bcsr        = csqnbar - csqn-cacsr-dcacsr+kmcsqn;
	ccsr        = kmcsqn*(csqn+cacsr+dcacsr);
	cacsr       = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;
	
	dcansr	    = dt*(qup1+qup2-qtr1*vjsr/vnsr-qtr2*vcsr/vnsr);
 	cansr	   += dcansr;
 	
	dnasss	    = dt*((-(3*inacass)*acap)/((vsss)*zna*frdy)-qdiffna); 
	nasss      += dnasss;
	
	dnassl	    = dt*((-(3*inak+ina+inal+3*inaca+ifna+inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
	nassl	   += dnassl;
	
	dnai        = dt*(qgapna*vssl/vmyo);
	nai        += dnai;
	
	dki	        = dt*((-iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	ki         += dki;
	
	dcai		= dt*(-(qup2)*vnsr/vmyo+qgap*vssl/vmyo+(qrel2)*vcsr/vmyo);
	trpn        = trpnbar*(cai/(cai+kmtrpn));
	cmdn		= cmdnbar*(cai/(cai+kmcmdn));
	catotal		= trpn+cmdn+dcai+cai;
	bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cai		    = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;  
	
	caavg       = (casss*vsss+cassl*vssl+cai*vmyo)/(vsss+vmyo+vssl);
	
 	camkbound	= camk0*(1-camktrap)*1/(1+(kmcam/casss));
	camktrap	= dt*(alphacamk*camkbound*(camkbound+camktrap)-betacamk*camktrap) + camktrap;
	camkactive	= camkbound+camktrap; 
	
}             

void printtofile ()
{
	count    += 1;
	
	//PRINT LAST 5 BEATS
	if (count>=10 && t>=(BCL*(beats-5)))			
	{
		count=0;
		fprintf(ap,"%f\t%f\t%f\n", t-BCL*(beats-5), v, caavg);
   	}	
	
}
