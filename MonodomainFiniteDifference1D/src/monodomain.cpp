#include "../include/monodomain.h"

struct monodomain_solver* new_monodomain_solver ()
{
	struct monodomain_solver *solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
	
	solver->sv = NULL;
	solver->stim_current = NULL;
	solver->vm = NULL;

	return solver; 
}

void free_monodomain_solver (struct monodomain_solver *solver)
{
	if (solver->sv)
		delete [] solver->sv;

	if (solver->stim_current)
		delete [] solver->stim_current;

	if (solver->vm)
		delete [] solver->vm;

	free(solver);
}

void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options)
{
	solver->dx = options->dx;
	solver->dt = options->dt;
	solver->tmax = options->tmax;
	solver->lmax = options->lmax;
	solver->Ncell = nearbyint(solver->lmax / solver->dx);
	solver->Niter = nearbyint(solver->tmax / solver->dt);

	solver->use_steady_state = options->use_steady_state;
	solver->sst_filename = options->sst_filename;

	// Allocate memory
	solver->sv = new double[solver->Ncell*Nodes]();
	solver->stim_current = new double[solver->Ncell]();
	solver->vm = new double[solver->Ncell]();

	// DEBUG
	//print_monodomain_solver(solver);
}

void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter)
{
	// Get the reference to the structures variables 
	double *sv = solver->sv;
	double *stim_current = solver->stim_current;
	double *vm = solver->vm;
	double dx = solver->dx;
	double dt = solver->dt;
	int Ncell = solver->Ncell;
	int Niter = solver->Niter;

	int print_rate = plotter->print_rate;
	int sst_rate = plotter->sst_rate;
	std::ofstream *plot_files = plotter->plot_files;
	int *plot_cell_ids = plotter->plot_cell_ids;

	#ifdef OUTPUT_VTK
	printf("[Monodomain] Saving output files in VTK format\n");
	#endif

	std::string sst_filename = solver->sst_filename; 
	
	// Initial conditions configuration
	if (solver->use_steady_state)
		read_initial_conditions_from_file(sv,Ncell,Nodes,sst_filename);
	else
		compute_initial_conditions(sv,Ncell,Nodes);

	double start, finish, elapsed;
	GET_TIME(start);	

	for (int k = 0; k <= Niter; k++)
	{
		double t = dt*k;

		print_progress(k,Niter);	

		if (k % print_rate == 0)
		{
			#ifdef OUTPUT_VTK
			write_VTK_to_file(sv,dx,Ncell,Nodes,k);
			#endif
			write_plot_data(plot_files,t,sv,Ncell,Nodes,plot_cell_ids);
		}

		if (k % (sst_rate-1) == 0)
		{
			write_steady_state_to_file(sv,Ncell,Nodes);
		}

		compute_stimulus(stim,stim_current,t,Ncell,dx);
		//print_stimulus(stim_current,Ncell,dx);
		
		solve_diffusion(sv,dx,dt,Ncell,Nodes,vm);

		update_state_vector(sv,vm,Ncell,Nodes);		
		
		solve_reaction(sv,stim_current,t,Ncell,Nodes,dt);

	}

	GET_TIME(finish);
	elapsed = finish - start;
	printf("%s\n",PRINT_LINE);
	printf("Elapsed time = %.10lf\n",elapsed);
	printf("%s\n",PRINT_LINE);
}

void solve_diffusion (double *sv, const double dx, const double dt,\
 			const int np, const int nodes,\
			double *vm)
{
	static const double D = 2.5e-04; // Paper: D = 2.5e-04
	static const double Cm = 12.0;
	static const double r = (dt*D)/(dx*dx);
	
	//#pragma omp parallel 
	for (int i = 0; i < np; i++)
	{
		// Boundary node
		if (i == 0)
			vm[i] = (1.0 - r)*sv[i*nodes] + (r)*sv[(i+1)*nodes];
		// Boundary node		
		else if (i == np-1)
			vm[i] = (1.0 - r)*sv[i*nodes] + (r)*sv[(i-1)*nodes];
		// Interior node
		else
			vm[i] = (1.0 - 2.0*r)*sv[i*nodes] + (r)*sv[(i+1)*nodes] + (r)*sv[(i-1)*nodes];
	} 
	
}

void update_state_vector (double *sv, const double *vm,\
			  const int np, const int nodes)
{
	#pragma omp parallel for
	for (int i = 0; i < np; i++)
	{
		sv[i*nodes] = vm[i];
	}
}

void solve_reaction (double *sv, double *stims, const double t,\
		     const int np, const int nodes, const double dt)
{
	#pragma omp parallel for
	for (int i = 0; i < np; i++)
	{
		double V_old = sv[i*nodes];
		double m_old = sv[i*nodes+1];
		double h_old = sv[i*nodes+2];
		double n_old = sv[i*nodes+3];

		sv[i*nodes] = V_old + dt*dvdt(V_old,m_old,h_old,n_old,stims[i]);
		sv[i*nodes+1] = m_old + dt*dmdt(V_old,m_old);
		sv[i*nodes+2] = h_old + dt*dhdt(V_old,h_old);
		sv[i*nodes+3] = n_old + dt*dndt(V_old,n_old);
		
	}
}

void print_monodomain_solver (const struct monodomain_solver *solver)
{

	printf("[Monodomain] dx = %.10lf cm\n",solver->dx);
	printf("[Monodomain] dt = %.10lf ms\n",solver->dt);
	printf("[Monodomain] tmax = %.10lf ms\n",solver->tmax);
	printf("[Monodomain] lmax = %.10lf cm\n",solver->lmax);
	printf("[Monodomain] Number of cells = %d\n",solver->Ncell);
	printf("[Monodomain] Number of iterations = %d\n",solver->Niter);
	printf("[Monodomain] Number of ODE equations = %d\n",Nodes);
}
