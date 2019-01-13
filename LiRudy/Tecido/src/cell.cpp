#include "../include/cell.h"

// GEOMETRY
double ageo,acap,vcell,vmyo,vnsr,vjsr,vsss,vcsr,vssl,vmito;

double __dt;
double __tmax;

Cell::Cell ()
{
    
}

void Cell::timestep (double dt)
{
    dvdtclock += dt;
}

void setTimestep (double dt, double tmax)
{
    __dt = dt;
    __tmax = tmax;
}

// Compute the geometry related variables 
void computeGeometrics ()
{
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
}

void Cell::setCell ()
{
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
}


void Cell::comp_conc ()
{
	qdiff       = (casss-cassl)/sstau;  
	qgap        = (cassl-cai)/gaptau;  
    qdiffna     = (nasss-nassl)/sstau;
    qgapna      = (nassl-nai)/gaptau;

    //printf("qdiff = %lf\n",qdiff);
	//printf("qgap = %lf\n",qgap);
	//printf("qdiffna = %lf\n",qdiffna);
	//printf("qgapna = %lf\n",qgapna);
    
	dcasss		= __dt*(-(ical-2*inacass)*acap/(vsss*2.0*frdy)+(qrel1+qip3)*vjsr/vsss-qdiff);
	bsss        = 1/(1+(bsrbar*kmbsr/pow(kmbsr+casss,2))+(bslbar*kmbsl/pow(kmbsl+casss,2)));
	casss      += bsss*dcasss;
    //printf("casss = %lf\n",casss);
	
	dcassl		= __dt*(-(qup1)*vnsr/vssl+qdiff*vsss/vssl-qgap-(icat+ipca+icab-2*inaca)*acap/(vssl*2.0*frdy));
	trpn        = trpnbar1*(cassl/(cassl+kmtrpn));
	cmdn		= cmdnbar1*(cassl/(cassl+kmcmdn));
	catotal		= trpn+cmdn+dcassl+cassl;
	bmyo		= cmdnbar1+trpnbar1-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar1*kmcmdn)+cmdnbar1*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cassl		= (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;   
 	//printf("casss = %lf\n",cassl);

	dcajsr		= __dt*(qtr1-qrel1-qip3);
	csqn1       = csqnbar1*(cajsr/(cajsr+kmcsqn));
	bjsr        = csqnbar1 - csqn1-cajsr-dcajsr+kmcsqn;
	cjsr        = kmcsqn*(csqn1+cajsr+dcajsr);
	cajsr       = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
	//printf("casss = %lf\n",cajsr);

	dcacsr		= __dt*(qtr2-qrel2);
	csqn        = csqnbar*(cacsr/(cacsr+kmcsqn));
	bcsr        = csqnbar - csqn-cacsr-dcacsr+kmcsqn;
	ccsr        = kmcsqn*(csqn+cacsr+dcacsr);
	cacsr    = (sqrt(bcsr*bcsr+4*ccsr)-bcsr)/2;
	//printf("casss = %lf\n",cacsr);

	dcansr	    = __dt*(qup1+qup2-qtr1*vjsr/vnsr-qtr2*vcsr/vnsr);
 	cansr	   += dcansr;
    //printf("cansr = %lf\n",cansr);
 	
	dnasss	    = __dt*((-(3*inacass)*acap)/((vsss)*zna*frdy)-qdiffna); 
	nasss      += dnasss;
    //printf("nasss = %lf\n",nasss);
	
	dnassl	    = __dt*((-(3*inak+ina+inal+3*inaca+ifna+inab)*acap)/((vssl)*zna*frdy)+qdiffna*vsss/vssl-qgapna);
	nassl	   += dnassl;
    //printf("nassl = %lf\n",nassl);
	
	dnai        = __dt*(qgapna*vssl/vmyo);
	nai        += dnai;
    //printf("nai = %lf\n",nai);
	
	dki	        = __dt*((-iktot*acap)/((vmyo+vssl+vsss)*zk*frdy));
	ki         += dki;
    //printf("ki = %lf\n",ki);
	
	dcai		= __dt*(-(qup2)*vnsr/vmyo+qgap*vssl/vmyo+(qrel2)*vcsr/vmyo);
	trpn        = trpnbar*(cai/(cai+kmtrpn));
	cmdn		= cmdnbar*(cai/(cai+kmcmdn));
	catotal		= trpn+cmdn+dcai+cai;
	bmyo		= cmdnbar+trpnbar-catotal+kmtrpn+kmcmdn;
	cmyo		= kmcmdn*kmtrpn-catotal*(kmtrpn+kmcmdn)+(trpnbar*kmcmdn)+cmdnbar*kmtrpn;
	dmyo		= -kmtrpn*kmcmdn*catotal;
	cai		    = (2.0/3.0)*sqrt(bmyo*bmyo-3.0*cmyo)*cos(acos((9.0*bmyo*cmyo-2*bmyo*bmyo*bmyo-27*dmyo)/(2.0*pow((bmyo*bmyo-3.0*cmyo),1.5)))/3.0)-bmyo/3.0;  
	
	caavg       = (casss*vsss+cassl*vssl+cai*vmyo)/(vsss+vmyo+vssl);
	//printf("caavg = %lf\n",caavg);

 	camkbound	= camk0*(1-camktrap)*1/(1+(kmcam/casss));
	camktrap	= __dt*(alphacamk*camkbound*(camkbound+camktrap)-betacamk*camktrap) + camktrap;
	camkactive	= camkbound+camktrap; 
	//printf("camkactive = %lf\n",camkactive);

}         

void Cell::comp_qtr2 ()
{
	qtr2		= (cansr-cacsr)/tautr2;
}

void Cell::comp_qtr1 ()
{
	qtr1		= (cansr-cajsr)/tautr1;
}

void Cell::comp_qup2 ()
{
    dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	qup2		= 0.0026*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cai,1))-0.0042*cansr/nsrbar;
}

void Cell::comp_qup1 ()
{
    dkmplb		= dkmplbbar*camkactive/(kmcamk+camkactive);
	dqupcamk	= dqupcamkbar*camkactive/(kmcamk+camkactive); 
	qup1		= 0.0002*(dqupcamk+1)/(1+pow((kmup-dkmplb)/cassl,1))-0.00105*cansr/nsrbar;
}

void Cell::comp_qrel2 ()
{
	qgap  = (cassl-cai)/gaptau;  
    REL  = (-qup2*vnsr/vmyo + qgap*vssl/vmyo+ (qrel2)*vcsr/vmyo);    
    ireltau = 6*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cacsr));
    if (REL > 0)
        irelss  = 91*(1+1*(1/(1+pow((0.28/camkactive),8))))*(REL)/(1 + pow((1/cacsr),8));
    else 
        irelss = 0;
    qrel2 += __dt*((irelss-qrel2)/ireltau);
}

void Cell::comp_qrel1 ()
{
	qdiff  = (casss-cassl)/sstau;  
    REL  = -((ical)*acap/(vsss*2.0*frdy) - (qrel1 + qip3)*vjsr/vsss + qdiff);     
    ireltau = 2*(1+1*(1/(1+pow((0.28/camkactive),8))))/(1+(0.0123/cajsr));
    if (REL > 0)
        irelss  = 15*(1+1*(1/(1+pow((0.28/camkactive),8))))*REL/(1 + pow((1.0/cajsr),8));
    else 
        irelss = 0;
    qrel1 += __dt*((irelss-qrel1)/ireltau);
}

void Cell::comp_ip3 ()
{
    u += __dt*(casss*k2*(1-u) - k2a*u);
    POip3 = tauip3r*IP3*casss*(1-u)/((1+IP3*k0/k0a)*(1+casss*k1/k1a));
    qip3 = 10.920*(cajsr-casss)*(POip3);
}

void Cell::comp_itot ()
{
	if (stimtime >= 0.0 && stimtime < stimdur)
	{
		//printf("stimtime = %lf || stimdur = %lf\n",stimtime,stimdur);
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

void Cell::comp_istim (double t) 
{
	stimtime += __dt;
	if (t >= tstim)
	{
		stimtime = 0.0;					
		stimcount += 1;					
		if (stimcount < beats-1)  tstim += BCL;		
		else if (stimcount == beats-1) tstim += S2;	
		else tstim = __tmax + 1;				
		//if (stimcount < beats) printf ("S1 Beat %d at time = %.2f ms !\n", stimcount+1, t);
		//else if (stimcount == beats) printf ("S2 Beat at time = %.2f ms !\n", t);
	}
}

void Cell::comp_if ()
{
	yss       = 1/(1+exp((v+87)/9.5));
	ytau      = 2000/(exp(-(v+132)/10) + exp((v+57)/60));
	y      = yss - (yss-y)*exp(-__dt/ytau);
	ifna	  = 0.012*y*y*(v-ena);
	ifk       = 0.024*y*y*(v-ek);
	iftotal   = ifna + ifk;
}

void Cell::comp_ipca ()
{
    ipca = ipcabar/((kmpca/cassl)+1);
}

void Cell::comp_inak ()
{
    inak	= ibarnak*(1/(1+exp(-1*(v+92)*frdy/(R*temp))))*pow((nassl/(nassl+2.6)),3)*(ko/(ko+0.8));
}

void Cell::comp_inaca ()
{
	allo		= 1/(1+pow((kmcaact/(1.5*casss)),2));
	num		    = inacamax*(pow(nasss,3)*cao*exp(nu*v*frdy/(R*temp))-pow(nao,3)*1.5*casss*exp((nu-1)*v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v*frdy/(R*temp));
	denomterm1	= kmcao*pow(nasss,3)+pow(kmnao,3)*1.5*casss+pow(kmnai1,3)*cao*(1+1.5*casss/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nasss/kmnai1,3))+pow(nasss,3)*cao+pow(nao,3)*1.5*casss;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inacass  = 0.2*allo*deltaE;
	
	allo		= 1/(1+pow((kmcaact/(1.5*cassl)),2));
	num		    = inacamax*(pow(nassl,3)*cao*exp(nu*v*frdy/(R*temp))-pow(nao,3)*1.5*cassl*exp((nu-1)*v*frdy/(R*temp)));
	denommult	= 1+ksat*exp((nu-1)*v*frdy/(R*temp));
	denomterm1	= kmcao*pow(nassl,3)+pow(kmnao,3)*1.5*cassl+pow(kmnai1,3)*cao*(1+1.5*cassl/kmcai);
	denomterm2	= kmcai*pow(nao,3)*(1+pow(nassl/kmnai1,3))+pow(nassl,3)*cao+pow(nao,3)*1.5*cassl;
	deltaE		= num/(denommult*(denomterm1+denomterm2));
	inaca    = 0.8*allo*deltaE;
}

void Cell::comp_ik1 ()
{	
    k1ss      = 1/(1+exp((v+103-(2.9+ko*2.175))/10.15));
	gk1	      = 0.12*sqrt(ko);
	ik1	      = gk1*k1ss*(v-ek);
}

void Cell::comp_iks ()
{
	eks	    = (R*temp/frdy)*log((ko+prnak*nao)/(ki+prnak*nassl));
	gks	    = 0.053*(1+0.6/(1+pow((0.000038/cassl),1.4)));
	xsss	= 1/(1+exp(-(v-9)/13.7));
	xs1tau	= 200/(exp(-(v+10)/6) + exp((v-62)/55));
	xs2tau	= 1500+ 350/(exp(-(v+10)/4) + exp((v-90)/58));
	xs1	= xsss-(xsss-xs1)*exp(-__dt/xs1tau);
	xs2	= xsss-(xsss-xs2)*exp(-__dt/xs2tau);
	iks	= gks*xs1*xs2*(v-eks);
}

void Cell::comp_ikr ()
{
	gkr	    = 0.0326*sqrt(ko/5.4);
	xrss	= 1/(1+exp(-(v)/15));
	xrtau   = 400.0/(1.0+exp(v/10.0)) + 100.0;
	rkr	    = 1/(1+exp((v)/35));
	xr	    = xrss-(xrss-xr)*exp(-__dt/xrtau);
	ikr	    = gkr*xr*rkr*(v-ek);
}

void Cell::comp_ito1 ()
{
	atau	= 1/(25*exp((v-82)/18)/(1+exp((v-82)/18))+25*exp(-(v+52)/18)/(1+exp(-(v+52)/18)));
	itau	= 2.86+ 1/(exp(-(v+125)/15)*0.1 + 0.1*exp((v+2)/26.5));
	i2tau	= 21.5+ 1/(exp(-(v+138.2)/52)*0.005 + 0.003*exp((v+18)/12.5));
	ass	    = 1/(1+exp(-(v-8.9)/10.3));
	iss	    = 1/(1+exp((v+30)/11));
	i2ss	= iss;
	a	    = ass-(ass-a)*exp(-__dt/atau);
	i	    = iss-(iss-i)*exp(-__dt/itau);
	i2	    = i2ss-(i2ss-i2)*exp(-__dt/i2tau);
	itos    = gtos*a*i*i2*(v-ek);
	itof    = gtof*(v-ek)/(1+exp(-(v-3)/19.8));
	ito1	= itos + itof;
}

void Cell::comp_icab ()
{
	icab	= pcab*zca*zca*((v*frdy*frdy)/(R*temp))*((gacai*cassl*exp((zca*v*frdy)/(R*temp))-gacao*cao)/(exp((zca*v*frdy)/(R*temp))-1));
}

void Cell::comp_icat ()
{
	bss	    = 1/(1+ exp (-(v+30)/7));
	gss	    = 1/(1+exp((v+61)/5));
	taub	= 1/(1.068*exp((v+16.3)/30)+1.068*exp(-(v+16.3)/30));
	taug    = 1/(0.015*exp(-(v+71.7)/83.3)+0.015*exp((v+71.7)/15.4));
	b	= bss-(bss-b)*exp(-__dt/taub);
	g	= gss-(gss-g)*exp(-__dt/taug);
	icat	= gcat*b*g*(v-eca);
}

// Pode ser que precise colocar as corrente 'ical' como variavel da struct Cell 
void Cell::comp_ical ()
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
	d		= dss-(dss-d)*exp(-__dt/dtau);
	f		= fss-(fss-f)*exp(-__dt/ftau);
	f2		= f2ss-(f2ss-f2)*exp(-__dt/f2tau);
	fca		= fcass-(fcass-fca)*exp(-__dt/fcatau);
	fca2		= fca2ss-(fca2ss-fca2)*exp(-__dt/fca2tau);
	ical		= d*f*f2*fca*fca2*ibarca;	
}

void Cell::comp_inab ()
{
    inab    = pnab*frdy*((frdy*v)/(R*temp))*(nassl*exp((frdy*v)/(R*temp)) - nao)/(exp((frdy*v)/(R*temp))-1);     
}

void Cell::comp_inal ()
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
	ml	  = mlss-(mlss-ml)*exp(-__dt/mltau);
	ml3     = ml3ss-(ml3ss-ml3)*exp(-__dt/ml3tau);
	hl	  = hlss-(hlss-hl)*exp(-__dt/hltau);
	hl3     = hl3ss-(hl3ss-hl3)*exp(-__dt/hl3tau);
	jl	  = jlss-(jlss-jl)*exp(-__dt/jltau);
	jl3     = jl3ss-(jl3ss-jl3)*exp(-__dt/jl3tau);
	inal2   = gnal2*ml*hl*jl*(v-ena);
	inal3   = gnal3*ml3*hl3*jl3*(v-ena);
	inal    = inal2 + inal3; 
}

void Cell::comp_ina ()
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
	m	= mss-(mss-m)*exp(-__dt/mtau);
	h	= hss-(hss-h)*exp(-__dt/htau);
	j	= jss-(jss-j)*exp(-__dt/jtau);
	ina	= gna*pow(m,3)*h*j*(v-ena);
}

void Cell::comp_revs ()
{
    eca	= (R*temp/(zca*frdy))*log(cao/cassl);
	ena	= (R*temp/frdy)*log(nao/nassl);
	ek	= (R*temp/frdy)*log(ko/ki);
}