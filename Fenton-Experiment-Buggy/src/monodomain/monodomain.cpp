#include "monodomain.h"

static inline double ALPHA (double dt, double h, double D) 
{
    return ( D * dt) / ( pow(h,2.0) );
}

/*
static inline double ALPHA (const double beta, const double cm, const double d, const double dx, const double dt) 
{
    //return ((( M_PI * (beta) * (cm) * (d) * (d) * (dx) ) / ( (dt) * (4.0) ) ) * UM2_TO_CM2);		// Cylinder control volume
	return ((( (beta) * (cm) * (dx) * (dx) * (dx) ) / ( (dt) ) ) * UM2_TO_CM2);						// Cubic control volume	
}
*/

static inline double GAMMA (const double d, const double sigma_c, const double dx) 
{
    //return (M_PI * d * d * sigma_c) / (4.0 * dx);													// Cylinder control volume
	return ( (sigma_c) * (dx) * (dx) ) / (dx);														// Cubic control volume
}

struct monodomain_solver* new_monodomain_solver ()
{
	struct monodomain_solver *solver = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));
	
	solver->sv = NULL;
	solver->stim_current = NULL;
	solver->vm = NULL;
	solver->cell_mask = NULL;
	solver->activation_time = NULL;
	solver->max_dvdt = NULL;

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

	if (solver->cell_mask)
		delete [] solver->cell_mask;
	
	if (solver->max_dvdt)
		delete [] solver->max_dvdt;

	if (solver->activation_time)
		delete [] solver->activation_time;

	free(solver);
}

void configure_solver_from_options (struct monodomain_solver *solver, struct user_options *options)
{
	solver->num_threads = options->num_threads;
	solver->dx = options->dx;
	solver->dt = options->dt;
	solver->tmax = options->tmax;
	solver->lmax = options->lmax;
	solver->diameter = options->diameter;
	solver->num_cell_div = options->num_cell_div;
	solver->Ncell = nearbyint(solver->lmax / solver->dx);
	solver->Niter = nearbyint(solver->tmax / solver->dt);

	solver->sigma_c = options->sigma_c;
	solver->G_gap = options->G_gap;

	solver->use_steady_state = options->use_steady_state;
	solver->sst_filename = options->sst_filename;
	solver->calc_activation_time = options->calc_activation_time;

	// Allocate memory
	solver->sv = new double[solver->Ncell*Nodes]();
	solver->stim_current = new double[solver->Ncell]();
	solver->vm = new double[solver->Ncell]();
	solver->cell_mask = new bool[solver->Ncell];
	if (solver->calc_activation_time)
	{
		solver->activation_time = new double[solver->Ncell];
		solver->max_dvdt = new double[solver->Ncell];

		configure_activation_time(solver->activation_time,solver->max_dvdt,solver->Ncell);
	}

	configure_cell_mask(solver->cell_mask,solver->Ncell,solver->num_cell_div);

	// DEBUG
	//print_monodomain_solver(solver);
}

void configure_cell_mask (bool cell_mask[], const int ncell, const int num_cell_div)
{
	// Homogenous model (without Ggap)
	if (num_cell_div == 0)
		memset(cell_mask,false,ncell);
	// Heterogenous model (with Ggap)
	else
	{
		for (int i = 0; i < ncell; i++)
		{
			if (i % num_cell_div == 0 && i != 0)
				cell_mask[i] = true;
			else
				cell_mask[i] = false;
		}
	}

}

void configure_activation_time(double *activation_time, double *max_dvdt, const int ncell)
{
	memset(activation_time,-1.0,sizeof(double)*ncell);
	memset(max_dvdt,__DBL_MIN__,sizeof(double)*ncell);
}

void solve_monodomain (struct monodomain_solver *solver,\
			struct stim_config *stim,\
			struct plot_config *plotter)
{
	// Constants
	const double Cm = 1.0;
	const double beta = 0.14;

	// Get the reference to the structures variables 
	double *sv = solver->sv;
	double *stim_current = solver->stim_current;
	double *vm = solver->vm;
	bool *mask = solver->cell_mask;
	double *at = solver->activation_time;
	double *max_dvdt = solver->max_dvdt;
	double dx = solver->dx;
	double diameter = solver->diameter;
	double dt = solver->dt;
	bool calc_activation_time = solver->calc_activation_time;
	
	int Ncell = solver->Ncell;
	int Niter = solver->Niter;

	double sigma_c = solver->sigma_c;
	double G_gap = solver->G_gap;

	//double alpha = ALPHA(beta,Cm,diameter,dx,dt);
	static const double D = 2.5e-04;
	double alpha = ALPHA(D,dx,dt);
	double gamma = GAMMA(diameter,sigma_c,dx);

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
		
		solve_diffusion(sv,vm,mask,alpha,gamma,G_gap,Ncell,Nodes);

		update_state_vector(sv,vm,Ncell,Nodes);		
		
		solve_reaction(sv,stim_current,t,Ncell,Nodes,dt);

		if (calc_activation_time)
			calculate_derivative(at,max_dvdt,sv,vm,Ncell,Nodes,dt,t);


	}

	if (calc_activation_time)
	{
		calculate_propagation_velocity(at,plot_cell_ids,dx,Ncell);
		write_activation_times(at,max_dvdt,Ncell);
	}

	GET_TIME(finish);
	elapsed = finish - start;
	printf("%s\n",PRINT_LINE);
	printf("Elapsed time = %.10lf\n",elapsed);
	printf("%s\n",PRINT_LINE);
}



void solve_diffusion (const double *sv, double *vm, const bool *mask, const double alpha, const double gamma, const double G_gap, const int ncell, const int nodes)
{

	// Explicit:
	/*
	#pragma omp parallel for
	for (int i = 0; i < ncell; i++)
	{
		double west_flux, east_flux;

		// Case 1: First volume
		if (i == 0)
		{
			double multiplier = calculate_flux(mask[i+1],gamma,G_gap);
			east_flux = (-multiplier) * (sv[nodes*(i+1)] - sv[nodes*i]);
			west_flux = 0.0;
		}
		// Case 2: Last volume
		else if (i == ncell-1)
		{
			double multiplier = calculate_flux(mask[i-1],gamma,G_gap);
			east_flux = 0.0;
			west_flux = (-multiplier) * (sv[nodes*(i)] - sv[nodes*(i-1)]);
		}
		// Case 3: Middle volume
		else
		{
			double multiplier_east = calculate_flux(mask[i+1],gamma,G_gap);
			double multiplier_west = calculate_flux(mask[i-1],gamma,G_gap);
			east_flux = (-multiplier_east) * (sv[nodes*(i+1)] - sv[nodes*i]);
			west_flux = (-multiplier_west) * (sv[nodes*(i)] - sv[nodes*(i-1)]);
		}
		
		vm[i] = sv[nodes*i] + ( (west_flux - east_flux) / alpha );
		
	}	
	*/
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

double calculate_flux (const bool type, const double gamma, const double G_gap)
{
	if (type == true)
		return G_gap;
	else
		return gamma;
}

void calculate_derivative (double *at, double *max_dvdt, double *sv, double *vms, const int ncell, const int nodes, const double dt, const double cur_time)
{
	//#pragma omp parallel for
	for (int i = 0; i < ncell; i++)
	{
		double v_new = sv[i*nodes];
		double v_old = vms[i];
		double dvdt = (v_new - v_old) / dt;

		if (dvdt > max_dvdt[i])
		{
			max_dvdt[i] = dvdt;
			at[i] = cur_time;
		}
	}

}

void calculate_propagation_velocity (const double *at, const int *plot_ids, const double dx, const int ncell)
{
	const int offset = 5;
	double plot_velocity[5];

	for (int i = 0; i < 5; i++)
	{
		int cell_index = plot_ids[i];

		double velocity = (offset*dx) / (at[cell_index+offset] - at[cell_index]);

		plot_velocity[i] = velocity * UM_PER_MS_TO_M_PER_S;
	}

	write_propagation_velocity(plot_velocity,plot_ids);
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

	printf("[Monodomain] num_threads = %d\n",solver->num_threads);
	printf("[Monodomain] dx = %.10lf cm\n",solver->dx);
	printf("[Monodomain] dt = %.10lf ms\n",solver->dt);
	printf("[Monodomain] tmax = %.10lf ms\n",solver->tmax);
	printf("[Monodomain] lmax = %.10lf cm\n",solver->lmax);
	printf("[Monodomain] sigma_c = %.10lf mS/cm\n",solver->sigma_c);
	printf("[Monodomain] G_gap = %.10lf uS\n",solver->G_gap);
	printf("[Monodomain] Number of cells = %d\n",solver->Ncell);
	printf("[Monodomain] Number of cell divisions = %d\n",solver->num_cell_div);
	printf("[Monodomain] Number of iterations = %d\n",solver->Niter);
	printf("[Monodomain] Number of ODE equations = %d\n",Nodes);

	//exit(EXIT_SUCCESS);
}
