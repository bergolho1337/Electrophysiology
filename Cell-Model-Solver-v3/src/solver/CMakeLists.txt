set(SOURCE_FILES_SOLVER
        output_utils.h output_utils.c
        config/config_common.c config/config_common.h
        config/stim_config_hash.c config/stim_config_hash.h
        config/stim_config.h config/stim_config.c
        config/config_parser.c config/config_parser.h
        ode_solver.c ode_solver.h
        constants.h)

ADD_LIBRARY(solvers STATIC ${SOURCE_FILES_SOLVER})
