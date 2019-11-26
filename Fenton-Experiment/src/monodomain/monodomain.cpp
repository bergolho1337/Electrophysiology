#include "monodomain.h"

static inline double ALPHA (const double beta, const double cm, const double d, const double dx, const double dt) 
{
    return (beta * cm * M_PI * d * d * dx) / (4.0 * dt);
}

struct monodomain_solver* new_monodomain_solver ()
{
	struct monodomain_solver *solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
	
	solver->sv = NULL;
	solver->stim_current = NULL;
	solver->vm = NULL;
	solver->link_mask = NULL;
	
	solver->the_matrix = new_tridiagonal_matrix_thomas();

	return solver; 
}

void free_monodomain_solver (struct monodomain_solver *solver)
{
	if (solver->sv)
		delete [] solver->sv;

	if (solver->stim_current)
		delete [] solver->stim_current;

	if (solver->link_mask)
		delete [] solver->link_mask;

	if (solver->vm)
		delete [] solver->vm;
	
	if (solver->the_matrix)
		free_tridiagonal_matrix_thomas(solver->the_matrix);

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
	solver->link_mask = new bool[solver->Ncell-1]();

	configure_cell_mask(solver);

	// DEBUG
	//print_monodomain_solver(solver);
}

void configure_cell_mask (struct monodomain_solver *solver)
{
	int Ncell = solver->Ncell;
	double dx = solver->dx;
	int num_cell_divisions = nearbyint(CELL_LENGTH / dx);
	printf("[Monodomain] There will be %d divisions for each cell\n",num_cell_divisions);

	static const bool USE_HOMOGENOUS_CONDUCTIVITY = true;

	// Homogenous model
	if (USE_HOMOGENOUS_CONDUCTIVITY)
	{
		memset(solver->link_mask,false,sizeof(bool)*(Ncell-1));
		solver->use_homogenous_conductivity = true;
	}
	// Heterogenous model
	else
	{
		for (int i = 0; i < Ncell-1; i++)
		{
			if (i % num_cell_divisions == 0 && i != 0)
				solver->link_mask[i] = true;
			else
				solver->link_mask[i] = false;
		}
		solver->use_homogenous_conductivity = false;
	}

}

void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter)
{
	// Get the reference to the structures variables 
	double *sv = solver->sv;
	double *stim_current = solver->stim_current;
	double *vm = solver->vm;
	bool *link_mask = solver->link_mask;
	double dx = solver->dx;
	double dt = solver->dt;
	double lmax = solver->lmax;
	double bcl = stim->start_period;
	int Ncell = solver->Ncell;
	int Niter = solver->Niter;
	bool use_homogenous_conductivity = solver->use_homogenous_conductivity;

	// Setting parameters
	double d = 0.001665;										// {cm} DO NOT CHANGE THIS VALUE !!!!
	double sigma = 4.0;											// {mS/cm}
	double G = 0.00831878940731400028;							// {uS}
	//double sigma = 2.0;
	//double G = 0.00839287531806615951;
	double beta = 0.14;											// {cm^-1}
	double cm = 1.0;											// {uF/cm^2}
	double sigma_bar = calc_homogenous_conductivity(sigma,G);	// {mS/cm}
	double alpha = ALPHA(beta,cm,d,dx,dt);
	
	print_monodomain_solver(solver,d,sigma,G,sigma_bar,beta,cm,use_homogenous_conductivity);

	// CONVERTING UNITS: Units should be given in {S/cm} and {S} for the solver
	sigma *= 1.0e-03;
	sigma_bar *= 1.0e-03;
	G *= 1.0e-06;

	int print_rate = plotter->print_rate;
	int sst_rate = plotter->sst_rate;
	std::ofstream *plot_files = plotter->plot_files;
	int *plot_cell_ids = plotter->plot_cell_ids;

	std::string sst_filename = solver->sst_filename; 
	
	#ifdef OUTPUT_VTK
	printf("[Monodomain] Saving output files in VTK format\n");
	#endif

	// Assembly matrix
	assembly_discretization_matrix(solver,sigma,G,sigma_bar,alpha,beta,cm,d,dx,dt,use_homogenous_conductivity);

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

		if (k % (sst_rate-1) == 0 && !solver->use_steady_state)
		{
			write_steady_state_to_file(sv,bcl,lmax,Ncell,Nodes);
		}

		compute_stimulus(stim,stim_current,t,Ncell,dx);
		//print_stimulus(stim_current,Ncell,dx);
		
		assembly_load_vector(solver,alpha);

		//solve_diffusion(sv,vm,link_mask,beta,cm,sigma,G,dx,d,dt,Ncell,Nodes);
		solve_diffusion_v2(solver);

		update_state_vector(sv,vm,Ncell,Nodes);		
		
		solve_reaction(sv,stim_current,t,Ncell,Nodes,dt);

	}

	GET_TIME(finish);
	elapsed = finish - start;
	printf("%s\n",PRINT_LINE);
	printf("Elapsed time = %.10lf\n",elapsed);
	printf("%s\n",PRINT_LINE);
}

void solve_diffusion (const double *sv, double *vm, bool *cell_mask, const double beta, const double cm, const double sigma, const double G, const double dx, const double d, const double dt, const int ncell, const int nodes)
{
	// Explicit: Finite Volume Method
	#pragma omp parallel for
	for (int i = 0; i < ncell; i++)
	{
		double total_flux = 0.0;
		double west_flux = 0.0;
		double east_flux = 0.0;

		// Case 1: First volume
		if (i == 0)
		{
			if (cell_mask[i+1] == true)
				east_flux = -G * (sv[nodes*(i+1)] - sv[nodes*i]);
				//east_flux = -sigma * ( (sv[nodes*(i+1)] - sv[nodes*i]) / dx ) * (M_PI * d * d / 4.0);
			else
				east_flux = -sigma * ( (sv[nodes*(i+1)] - sv[nodes*i]) / dx ) * (M_PI * d * d / 4.0);
		}
		// Case 2: Last volume
		else if (i == ncell-1)
		{
			if (cell_mask[i-1] == true)
				west_flux = -G * (sv[nodes*i] - sv[nodes*(i-1)]);
				//west_flux = -sigma * ( (sv[nodes*i] - sv[nodes*(i-1)]) / dx ) * (M_PI * d * d / 4.0);
			else
				west_flux = -sigma * ( (sv[nodes*i] - sv[nodes*(i-1)]) / dx ) * (M_PI * d * d / 4.0);
		}
		// Case 3: Middle volume
		else
		{
			if (cell_mask[i+1] == true)
				east_flux = -G * (sv[nodes*(i+1)] - sv[nodes*i]);
				//east_flux = -sigma * ( (sv[nodes*(i+1)] - sv[nodes*i]) / dx ) * (M_PI * d * d / 4.0);
			else
				east_flux = -sigma * ( (sv[nodes*(i+1)] - sv[nodes*i]) / dx ) * (M_PI * d * d / 4.0);
			
			if (cell_mask[i-1] == true)
				west_flux = -G * (sv[nodes*i] - sv[nodes*(i-1)]);
				//west_flux = -sigma * ( (sv[nodes*i] - sv[nodes*(i-1)]) / dx ) * (M_PI * d * d / 4.0);
			else
				west_flux = -sigma * ( (sv[nodes*i] - sv[nodes*(i-1)]) / dx ) * (M_PI * d * d / 4.0);
		}

		total_flux = east_flux - west_flux;

		vm[i] = sv[nodes*i] + ( (-total_flux / (beta*cm*M_PI*d*d*dx/4.0) * dt) );
	}	

	// Implicit: Solve the linear system
	//x = sparseSolver.solve(b);
}

void solve_diffusion_v2 (struct monodomain_solver *solver)
{
	uint32_t ncell = solver->Ncell;
	double *a = solver->the_matrix->a;
	double *b = solver->the_matrix->b;
	double *c = solver->the_matrix->c;
	double *d = solver->the_matrix->d;
	double *c_star = solver->the_matrix->c_star;
	double *d_star = solver->the_matrix->d_star;
	double *vms = solver->vm;							// Output

	solve_linear_system_using_thomas(a,b,c,d,c_star,d_star,vms,ncell);
}

void solve_linear_system_using_thomas (const double *a, const double *b, const double *c, const double *d,\
			  double *c_star, double *d_star, double *vms, const uint32_t ncell)
{
	uint32_t n = ncell - 1;		// Index starts at 0 not 1 ...
	
	// First sweep
	c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for (uint32_t i = 1; i < n; i++) 
    {
        c_star[i] = c[i] / (b[i] - a[i]*c_star[i-1]);
        d_star[i] = (d[i] - a[i]*d_star[i-1]) / (b[i] - a[i]*c_star[i-1]);
    }

	// Reverse sweep
    vms[n] = (d[n] - a[n]*d_star[n-1]) / (b[n] - a[n]*c_star[n-1]);
    for (uint32_t i = n; i-- > 0;) 
    {
        vms[i] = d_star[i] - c_star[i]*vms[i+1];
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

	//for (int i = 0; i < np; i++)
	//	printf("Volume = %d -- Vm = %g\n",i,sv[i*nodes]);
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

		// Explicit Euler
		//sv[i*nodes] = V_old + dt*dvdt(V_old,m_old,h_old,n_old,stims[i]);
		//sv[i*nodes+1] = m_old + dt*dmdt(V_old,m_old);
		//sv[i*nodes+2] = h_old + dt*dhdt(V_old,h_old);
		//sv[i*nodes+3] = n_old + dt*dndt(V_old,n_old);
		
		// Rush-Larsen
		sv[i*nodes] = V_old + dt*dvdt(V_old,m_old,h_old,n_old,stims[i]);
		sv[i*nodes+1] = dmdt_RL(V_old,m_old,dt);
		sv[i*nodes+2] = dhdt_RL(V_old,h_old,dt);
		sv[i*nodes+3] = dndt_RL(V_old,n_old,dt);
	}
	
}

struct tridiagonal_matrix_thomas* new_tridiagonal_matrix_thomas ()
{
	struct tridiagonal_matrix_thomas *result = (struct tridiagonal_matrix_thomas*)malloc(sizeof(struct tridiagonal_matrix_thomas));

	result->n = 0;
	result->a = NULL;
	result->b = NULL;
	result->c = NULL;
	result->d = NULL;
	result->c_star = NULL;
	result->d_star = NULL;

	return result;
}

void free_tridiagonal_matrix_thomas (struct tridiagonal_matrix_thomas *matrix)
{
	if (matrix->a)
		delete [] matrix->a;

	if (matrix->b)
		delete [] matrix->b;

	if (matrix->c)
		delete [] matrix->c;
	
	if (matrix->d)
		delete [] matrix->d;

	if (matrix->c_star)
		delete [] matrix->c_star;

	if (matrix->d_star)
		delete [] matrix->d_star;	
	
	free(matrix);
}

void assembly_discretization_matrix (struct monodomain_solver *solver,\
			 const double sigma_c, const double G_gap, const double sigma_bar, const double alpha, const double beta, const double cm,\
			 const double d, const double dx, const double dt,\
			 const bool use_homogenous_conductivity)
{
	struct tridiagonal_matrix_thomas *matrix = solver->the_matrix;

	// Allocate memory
	matrix->n = solver->Ncell;
	matrix->a = new double[solver->Ncell]();
	matrix->b = new double[solver->Ncell]();
	matrix->c = new double[solver->Ncell]();
	matrix->d = new double[solver->Ncell]();
	matrix->c_star = new double[solver->Ncell]();
	matrix->d_star = new double[solver->Ncell]();

	// Auxiliary variables
	uint32_t ncell = solver->Ncell;
	bool *link_mask = solver->link_mask;
	double gap_function_flux = G_gap;
	double citoplasm_flux;
	
	// HOMOGENOUS MODEL
	if (use_homogenous_conductivity)
		citoplasm_flux = (sigma_bar * M_PI * d * d) / (4.0 * dx);
	// HETEROGENOUS MODEL
	else
		citoplasm_flux = (sigma_c * M_PI * d * d) / (4.0 * dx);

	// Initalize the diagonal elements
	for (uint32_t i = 0; i < ncell; i++)
		matrix->b[i] = alpha;

	// Pass through each cell link
	for (uint32_t i = 0; i < ncell-1; i++)
	{
		bool link_type = link_mask[i];

		uint32_t west_cell_index = i;
		uint32_t east_cell_index = i+1;

		double flux;
		if (link_type)
			flux = gap_function_flux;
		else
			flux = citoplasm_flux;
		
		// East cell flux
		matrix->c[west_cell_index] = -flux;
		matrix->b[west_cell_index] += flux;

		// West cell flux
		matrix->a[east_cell_index] = -flux;
		matrix->b[east_cell_index] += flux;

	}

	// DEBUG
	//for (uint32_t i = 0; i < ncell; i++)
	//	printf("%u -- %20g %20g %20g\n",i,matrix->a[i],matrix->b[i],matrix->c[i]);
}

void assembly_load_vector (struct monodomain_solver *solver, const double alpha)
{
	uint32_t ncell = solver->Ncell;
	double *sv = solver->sv;
	double *d = solver->the_matrix->d;

	for (uint32_t i = 0; i < ncell; i++)
		d[i] = sv[i*Nodes] * alpha;
}

// Input: {mS/cm}, {uS}
// Output: {mS/cm}
double calc_homogenous_conductivity (const double sigma, const double G)
{
	const double l = 0.01;							// {cm}
	const double A = 235.6 * 1.0e-08;				// {cm^2}
	const double G_star = G * 1.0e-03;				// {mS}
	const double sigma_star = sigma;		// {mS/cm}

	// {mS/cm}
	double sigma_bar = l / ( (l/sigma_star) + (A/G_star) ); 

	return sigma_bar;
}

void print_monodomain_solver (const struct monodomain_solver *solver,\
		const double d, const double sigma, const double G, const double sigma_bar,const double beta,const double cm, const bool use_homogenous_conductivity)
{

	printf("[Monodomain] dx = %.10lf cm\n",solver->dx);
	printf("[Monodomain] dt = %.10lf ms\n",solver->dt);
	printf("[Monodomain] tmax = %.10lf ms\n",solver->tmax);
	printf("[Monodomain] lmax = %.10lf cm\n",solver->lmax);
	printf("[Monodomain] Number of cells = %d\n",solver->Ncell);
	printf("[Monodomain] Number of iterations = %d\n",solver->Niter);
	printf("[Monodomain] Number of ODE equations = %d\n",Nodes);
	printf("[Monodomain] Diameter = %g cm\n",d);
	printf("[Monodomain] Citoplasm conductivity = %g mS/cm\n",sigma);
	printf("[Monodomain] Gap junction conductance = %g uS\n",G);
	printf("[Monodomain] Homogenous conductivity = %g mS/cm\n",sigma_bar);
	printf("[Monodomain] beta = %g cm^-1\n",beta);
	printf("[Monodomain] Membrane capacitance = %g uF/cm^2\n",cm);
	if (use_homogenous_conductivity)
		printf("[Monodomain] Using HOMOGENOUS model\n");
	else
		printf("[Monodomain] Using HETEROGENOUS model\n");
}
