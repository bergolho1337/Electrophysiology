#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import numpy as np

NUM_THREADS = 4
DT_PDE = 0.02
SIMULATION_TIME = 500.0

UPDATE_LIBRARY_PATH = "shared_libs/libdefault_update_monodomain.so"

PRINT_RATE = 50
SAVE_RESULT_MAIN_FUNCTION = "save_as_vtp_purkinje"
SAVE_RESULT_LIBRARY_PATH = "shared_libs/libdefault_save_mesh.so"

SIGMA_X = 0.00005336
SIGMA_Y = 0.00005336
SIGMA_Z = 0.00005336
ASSEMBLY_MATRIX_LIBRARY_PATH = "shared_libs/libpurkinje_matrix_assembly.so"

LINEAR_SYSTEM_LIBRARY_PATH = "shared_libs/libdefault_linear_system_solver.so"

START_DISCRETIZATION = 100.0
PURKINJE_LIBRARY_PATH = "shared_libs/libdefault_purkinje.so"

DT_ODE = 0.02
ODE_SOLVER_LIBRARY_PATH = "shared_libs/libten_tusscher_3_sensibility.so"

STIMULUS_LIBRARY_PATH = "shared_libs/libdefault_stimuli.so"

EXTRA_LIBRARY_PATH = "shared_libs/libdefault_extra_data.so"

def write_main_section (file):
    file.write("[main]\n")
    file.write("num_threads = %d\n" % (NUM_THREADS))
    file.write("dt_pde = %g\n" % (DT_PDE))
    file.write("simulation_time = %g\n" % (SIMULATION_TIME))
    file.write("calc_activation_time = true\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n\n")

def write_update_monodomain_section (file):
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    file.write("library_file = %s\n\n" % (UPDATE_LIBRARY_PATH))

def write_save_result_section (file,atpi,Ko,Ki,Vm_modifier,GNa_modifier,GCaL_modifier,INaCa_modifier,num_biff):
    output_dir_name = "outputs/biff_%d/lucas:bifurcation-atpi_%g-Ko_%g-Ki_%g-Vm:mod_%g-GNa:mod_%g-GCaL:mod_%g-INaCa:mod_%g" % (num_biff,atpi,Ko,Ki,Vm_modifier,GNa_modifier,GCaL_modifier,INaCa_modifier)
    file.write("[save_result]\n")
    file.write("print_rate = %d\n" % (PRINT_RATE))
    file.write("output_dir = %s\n" % (output_dir_name))
    file.write("main_function = save_as_vtp_purkinje\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V\n")
    file.write("binary = false\n")
    file.write("compress = false\n")
    file.write("library_file = %s\n\n" % (SAVE_RESULT_LIBRARY_PATH))

def write_assembly_matrix_section (file):
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_fvm\n")
    file.write("sigma_x = %g\n" % (SIGMA_X))
    file.write("sigma_y = %g\n" % (SIGMA_Y))
    file.write("sigma_z = %g\n" % (SIGMA_Z))
    file.write("library_file = %s\n" % (ASSEMBLY_MATRIX_LIBRARY_PATH))
    file.write("main_function = purkinje_fibers_assembly_matrix\n\n")

def write_linear_system_section (file):
    file.write("[linear_system_solver]\n")
    file.write("tolerance = 1e-16\n")
    file.write("use_preconditioner = yes\n")
    file.write("max_iterations = 200\n")
    file.write("library_file = %s\n" % (LINEAR_SYSTEM_LIBRARY_PATH))
    file.write("main_function = conjugate_gradient\n\n")

def write_purkinje_section (file,num_biff):
    file.write("[purkinje]\n")
    file.write("name = Bifurcation Purkinje\n")
    file.write("start_discretization = %g\n" % (START_DISCRETIZATION))
    file.write("library_file = %s\n" % (PURKINJE_LIBRARY_PATH))
    file.write("main_function = initialize_purkinje_with_custom_mesh\n")
    file.write("network_file = networks/network_biff_%d.vtk\n\n" % (num_biff))

def write_ode_solver_section (file):
    file.write("[ode_solver]\n")
    file.write("dt_ode = %g\n" % (DT_ODE))
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = %s\n\n" % (ODE_SOLVER_LIBRARY_PATH))

def write_stimulus_section (file):
    file.write("[stim_purkinje_plain]\n")
    file.write("start = 1.0\n")
    file.write("duration = 2.0\n")
    file.write("current = -38.0\n")
    file.write("x_limit = 500.0\n")
    file.write("main_function = stim_if_x_less_than\n")
    file.write("library_file = %s\n\n" % (STIMULUS_LIBRARY_PATH))

def write_extra_data_section (file,atpi,Ko,Ki,Vm_modifier,GNa_modifier,GCaL_modifier,INaCa_modifier,Vm,M,H,J,Xr1,Xs,S,F,F2,D_inf,R_inf,Xr2_inf):
    file.write("[extra_data]\n")
    file.write("atpi = %g\n" % (atpi))
    file.write("Ko = %g\n" % (Ko))
    file.write("Ki = %g\n" % (Ki))
    file.write("Vm_modifier = %g\n" % (Vm_modifier))
    file.write("GNa_multiplicator = %g\n" % (GNa_modifier))
    file.write("GCaL_multiplicator = %g\n" % (GCaL_modifier))
    file.write("INaCa_multiplicator = %g\n" % (INaCa_modifier))
    file.write("sv_0 = %g\n" % (Vm))
    file.write("sv_1 = %g\n" % (M))
    file.write("sv_2 = %g\n" % (H))
    file.write("sv_3 = %g\n" % (J))
    file.write("sv_4 = %g\n" % (Xr1))
    file.write("sv_5 = %g\n" % (Xs))
    file.write("sv_6 = %g\n" % (S))
    file.write("sv_7 = %g\n" % (F))
    file.write("sv_8 = %g\n" % (F2))
    file.write("sv_9 = %g\n" % (D_inf))
    file.write("sv_10 = %g\n" % (R_inf))
    file.write("sv_11 = %g\n" % (Xr2_inf))
    file.write("main_function = set_extra_data_sensibility\n")
    file.write("library_file = %s\n\n" % (EXTRA_LIBRARY_PATH))


def main():

    if len(sys.argv) != 1:
                print("-------------------------------------------------------------------------")
                print("Usage:> python %s " % sys.argv[0])
                print("-------------------------------------------------------------------------")
                return 1

    os.chdir("inputs")
    for filename in glob.glob("*"):
        file = open(filename,"r")

        atpi = float(file.readline())
        Ko = float(file.readline())
        Ki = float(file.readline())
        Vm_mod = float(file.readline())
        GNa_mod = float(file.readline())
        GCaL_mod = float(file.readline())
        INaCa_mod = float(file.readline())

        Vm = float(file.readline())
        M = float(file.readline())
        H = float(file.readline())
        J = float(file.readline())
        Xr1 = float(file.readline())
        Xs = float(file.readline())
        S = float(file.readline())
        F = float(file.readline())
        F2 = float(file.readline())
        D_inf = float(file.readline())
        R_inf = float(file.readline())
        Xr2_inf = float(file.readline())

        file.close()

        for num_biff in range (1,5):
            output_filename = "/home/berg/Github/Electrophysiology/Broodie-Bifurcation-Ini-Generator/outputs/biff_%d/lucas:bifurcation-atpi_%g-Ko_%g-Ki_%g-Vm:mod_%g-GNa:mod_%g-GCaL:mod_%g-INaCa:mod_%g.ini" % (num_biff,atpi,Ko,Ki,Vm_mod,GNa_mod,GCaL_mod,INaCa_mod)
            output_file = open(output_filename,"w")

            write_main_section(output_file)
            write_update_monodomain_section(output_file)
            write_save_result_section(output_file,atpi,Ko,Ki,Vm_mod,GNa_mod,GCaL_mod,INaCa_mod,num_biff)
            write_assembly_matrix_section(output_file)
            write_linear_system_section(output_file)
            write_purkinje_section(output_file,num_biff)
            write_ode_solver_section(output_file)
            write_stimulus_section(output_file)
            write_extra_data_section(output_file,atpi,Ko,Ki,Vm_mod,GNa_mod,GCaL_mod,INaCa_mod,Vm,M,H,J,Xr1,Xs,S,F,F2,D_inf,R_inf,Xr2_inf)

            output_file.close()



if __name__ == "__main__":
        main()
