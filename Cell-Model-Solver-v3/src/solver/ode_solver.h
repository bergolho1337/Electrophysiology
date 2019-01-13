//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include "../utils/stop_watch.h"
#include "config/stim_config_hash.h"
#include "config/config_parser.h"

#include "../models_library/model_common.h"

struct ode_solver 
{

    void *handle;

    real max_dt;
    real min_dt;
    real rel_tol;
    real abs_tol;

    uint32_t *cells_to_solve;

    bool gpu;
    int gpu_id;

    real *sv;
    void *edo_extra_data;
    size_t extra_data_size;
    struct cell_model_data model_data;

    size_t pitch;

    //User provided functions
    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;

};

struct ode_solver* new_ode_solver();
void free_ode_solver(struct ode_solver *solver);

void configure_ode_solver_from_options (struct ode_solver *solver, struct user_options *options);
void init_ode_solver_with_cell_model (struct ode_solver* solver);
void set_ode_initial_conditions (struct ode_solver *solver);

void solve_celular_model (struct ode_solver *solver, struct user_options *options);
void solve_odes (struct ode_solver *solver, double cur_time, int ode_step, struct stim_config_hash *stim_configs);
void print_solver_info (struct ode_solver *solver, struct user_options *options);

void print_result(const struct ode_solver *solver, const struct user_options *configs, int count, double cur_time, bool save_in_binary, FILE *output_file);
void print_cell (const struct ode_solver *solver, FILE *output_file, double cur_time, bool save_in_binary);

#endif //MONOALG3D_EDO_SOLVER_H
