#include "../include/lirudy.h"

// OUTPUT FILE
int count;

LiRudy::LiRudy (int argc, char *argv[])
{
	dt = atof(argv[1]);
	tmax = atof(argv[2]);
	nthreads = atoi(argv[3]);
	n = tmax / dt;

	computeGeometrics();
	allocMem();
	setInitCond();
	setStimCells();
	setTimeSettings();
	setTimestep(dt,tmax);
}

void LiRudy::setTimeSettings ()
{
	for (int i = 0; i < NCELLS; i++)
    {
        cells[i].dvdtclock = 1000;
        cells[i].stimtime = 1000;
        cells[i].stimcount 	= -1;
    }
}

void LiRudy::setInitCond ()
{
	#pragma omp parallel for num_threads(nthreads)
	for (int i = 0; i < (int)cells.size(); i++)
    	cells[i].setCell();
}

void LiRudy::setStimCells ()
{
	for (int i = 0; i < NSC; i++)
    {
        cells[i].beats = 10;
        cells[i].BCL = 500;
        cells[i].S2 = 500;
        cells[i].tstim = 0;
    }
    for (int i = NSC; i < NCELLS; i++)
    {
        cells[i].beats = 10;
        cells[i].BCL = tmax + 500;
        cells[i].S2 = 500;
        cells[i].tstim = cells[i].BCL;
    }
}

void LiRudy::allocMem ()
{
	cells.assign(NCELLS,Cell());
}

void LiRudy::solve ()
{
	FILE *out = fopen("cell.dat","w+");
    // Time loop
	double t = 0;
    while (t <= tmax)
    {
		for (int i = 0; i < NCELLS; i++)
        	cells[i].timestep(dt);
	
		#pragma omp parallel for num_threads(nthreads)
		for (int i = 0; i < NCELLS; i++)
		{
			cells[i].comp_revs();
			cells[i].comp_ina ();
            cells[i].comp_inal ();
            cells[i].comp_inab ();
            cells[i].comp_ical ();
            cells[i].comp_icat ();
            cells[i].comp_icab ();
            cells[i].comp_ito1 ();
            cells[i].comp_ikr ();
            cells[i].comp_iks ();
            cells[i].comp_ik1 ();
            cells[i].comp_inaca ();
            cells[i].comp_inak ();
            cells[i].comp_ipca ();
            cells[i].comp_if ();
            cells[i].comp_istim (t);
            cells[i].comp_itot ();

            cells[i].comp_ip3 ();
            cells[i].comp_qrel1 ();
            cells[i].comp_qrel2 ();
            cells[i].comp_qup1 ();
            cells[i].comp_qup2 ();
            cells[i].comp_qtr1 ();
            cells[i].comp_qtr2 ();

            cells[i].comp_conc ();

			cells[i].dvdt	   = -cells[i].itot;
            cells[i].v	   	  += cells[i].dvdt*dt;
		}
		fprintf(out,"%.10lf\t%.10lf\n",t,cells[OUT_ID].v);
    	t += dt;
    }
    fclose(out);
}
