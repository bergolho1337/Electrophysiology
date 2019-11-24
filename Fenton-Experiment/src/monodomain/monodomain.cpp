#include "monodomain.h"

static inline double ALPHA (double dt, double h, double D) 
{
    return ( D * dt) / ( pow(h,2.0) );
}

struct monodomain_solver* new_monodomain_solver ()
{
	struct monodomain_solver *solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
	
	solver->sv = NULL;
	solver->stim_current = NULL;
	solver->vm = NULL;
	solver->cell_mask = NULL;

	return solver; 
}

void free_monodomain_solver (struct monodomain_solver *solver)
{
	if (solver->sv)
		delete [] solver->sv;

	if (solver->stim_current)
		delete [] solver->stim_current;

	if (solver->cell_mask)
		delete [] solver->cell_mask;

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
	solver->cell_mask = new bool[solver->Ncell]();

	configure_cell_mask(solver);

	// DEBUG
	//print_monodomain_solver(solver);
}

void configure_cell_mask (struct monodomain_solver *solver)
{
	int Ncell = solver->Ncell;
	double dx = solver->dx;
	int num_cell_divisions = nearbyint(CELL_LENGTH / dx);

	if (num_cell_divisions == 1)
	{
		memset(solver->cell_mask,false,sizeof(bool)*Ncell);
	}
	else if (num_cell_divisions > 1)
	{
		for (int i = 0; i < Ncell; i++)
		{
			if (i % num_cell_divisions == 0 && i != 0)
				solver->cell_mask[i] = true;
			else
				solver->cell_mask[i] = false;
		}
	}
	else
	{
		fprintf(stderr,"[-] ERROR! Invalid number for cell divisions! 'num_cell_divisions' = %d\n",num_cell_divisions);
		exit(EXIT_FAILURE);
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
	bool *cell_mask = solver->cell_mask;
	double dx = solver->dx;
	double dt = solver->dt;
	double lmax = solver->lmax;
	double bcl = stim->start_period;
	int Ncell = solver->Ncell;
	int Niter = solver->Niter;

	static const double d = 0.01;
	static const double sigma = 0.000035;
	static const double G = 0.000000628;
	static const double beta = 0.14;
	static const double cm = 1.0;
	double D = sigma / (beta*cm);
	double alpha = ALPHA(D,dx,dt);

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

		if (k % (sst_rate-1) == 0 && !solver->use_steady_state)
		{
			write_steady_state_to_file(sv,bcl,lmax,Ncell,Nodes);
		}

		compute_stimulus(stim,stim_current,t,Ncell,dx);
		//print_stimulus(stim_current,Ncell,dx);
		
		solve_diffusion(sv,vm,cell_mask,beta,cm,sigma,G,dx,d,dt,Ncell,Nodes);

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

/*
void assembly_matrix (Eigen::SparseMatrix<double> &A, const double h, const double dt,\
						const int ncell)
{
	printf("[Monodomain] Assembling linear system matrix\n");

	// Paper: D = ( sigma / (beta * cm) ) = 2.5e-04
	static const double D = 2.5e-06;

	// Non-zero coefficients
    std::vector< Eigen::Triplet<double> > coeff;

	for (int i = 0; i < ncell; i++)
	{
		// Case 1: First volume
		if (i == 0)
		{
			// Diagonal
			double diagonal_value = ALPHA(dt,h) - D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i,diagonal_value));

			// East flux
			double value = D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i+1,value));
		}
		// Case 2: Last volume
		else if (i == ncell-1)
		{
			// Diagonal
			double diagonal_value = ALPHA(dt,h) - D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i,diagonal_value));

			// West flux
			double value = D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i-1,value));
		}
		// Case 3: Middle volume
		else
		{
			// Diagonal
			double diagonal_value = ALPHA(dt,h) - 2.0*D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i,diagonal_value));

			// East flux
			double value = D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i+1,value));

			// West flux
			value = D*h;
			coeff.push_back(Eigen::Triplet<double>(i,i-1,value));
		}
	}

	// DEBUG
	//for (unsigned int i = 0; i < coeff.size(); i++)
	//	printf("(%d,%d) = %g\n",coeff[i].row(),coeff[i].col(),coeff[i].value());

	A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();
}
*/

/*
void assembly_load_vector (Eigen::VectorXd &b, const double *sv,
						const double h, const double dt, const int ncell, const int nodes)
{
	static const double D = 2.5e-04; 

	for (int i = 0; i < ncell; i++)
		b(i) = ALPHA(dt,h)*sv[nodes*i];
}
*/

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
