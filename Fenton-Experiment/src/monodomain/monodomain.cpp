#include "monodomain.h"

static inline double ALPHA (double dt, double h, double D) 
{
    return ( D * dt) / ( pow(h,2.0) * UM2_TO_CM2 );
}

static inline double CFL (double dt, double h, double D) 
{
    return ( D * dt / pow(h,2.0) );
}

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

	static const double D = 5.0e-03;
	//static const double D = 4.0e-04;
	double alpha = ALPHA(D,dx,dt);
	double cfl = CFL(D,dx,dt);
	if (cfl >= 0.5)
	{
		fprintf(stderr,"[Monodomain] ERROR! The CFL condition is not attended!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		printf("[Monodomain] CFL = %g\n",cfl);
	}

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
		
		solve_diffusion(sv,vm,alpha,Ncell,Nodes);

		update_state_vector(sv,vm,Ncell,Nodes);		
		
		solve_reaction(sv,stim_current,t,Ncell,Nodes,dt);

	}

	GET_TIME(finish);
	elapsed = finish - start;
	printf("%s\n",PRINT_LINE);
	printf("Elapsed time = %.10lf\n",elapsed);
	printf("%s\n",PRINT_LINE);
}

void solve_diffusion (const double *sv, double *vm, const double alpha, const int ncell, const int nodes)
{
	// Explicit:
	#pragma omp parallel for
	for (int i = 0; i < ncell; i++)
	{
		// Case 1: First volume
		if (i == 0)
			vm[i] = alpha*( sv[nodes*(i+1)] - sv[nodes*i] ) + sv[nodes*i];
		// Case 2: Last volume
		else if (i == ncell-1)
			vm[i] = alpha*( sv[nodes*(i-1)] - sv[nodes*i] ) + sv[nodes*i];
		// Case 3: Middle volume
		else
			vm[i] = alpha*( sv[nodes*(i-1)] + sv[nodes*(i+1)] - 2.0*sv[nodes*i] ) + sv[nodes*i];
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

		sv[i*nodes] = V_old + dt*dvdt(V_old,m_old,h_old,n_old,stims[i]);
		sv[i*nodes+1] = m_old + dt*dmdt(V_old,m_old);
		sv[i*nodes+2] = h_old + dt*dhdt(V_old,h_old);
		sv[i*nodes+3] = n_old + dt*dndt(V_old,n_old);
		
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
