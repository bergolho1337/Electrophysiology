#include "ini_parser/ini.h"
#include "solver/ode_solver.h"
#include "utils/logfile_utils.h"
#include "string/sds.h"
#include "solver/output_utils.h"

int main (int argc, char **argv) 
{
    if (argc-1 != 2)
    {
        display_usage(argv);
        exit(EXIT_FAILURE);
    }

    struct user_options *options;
    options = new_user_options ();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver ();

    // First we have to get the config file path
    get_config_file (argc, argv, options);

    if (options->config_file) 
    {
        // Here we parse the config file
        if (ini_parse (options->config_file, parse_config_file, options) < 0) 
        {
            fprintf (stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }

    // The command line options always overwrite the config file
    parse_options (argc, argv, options);

    // Create the output dir and the logfile
    if (options->out_dir_name) 
    {
        sds buffer = sdsnew ("");
        create_dir_if_no_exists (options->out_dir_name);
        buffer = sdscatfmt (buffer, "%s/outputlog.txt", options->out_dir_name);
        open_logfile (buffer);
        sdsfree (buffer);
    }

    configure_ode_solver_from_options (ode_solver, options);
    
    solve_celular_model(ode_solver,options);

    free_ode_solver (ode_solver);
    free_user_options (options);

    close_logfile ();

    return EXIT_SUCCESS;
}
