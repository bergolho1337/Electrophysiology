SET(LIB_TYPE MODULE) #set the lib type

############## MITCHELL SHAEFFER 2002 ##############################
set(MODEL_FILE_CPU "mitchell_shaeffer.c")
set(COMMON_HEADERS model_common.h mitchell_shaeffer.h)

ADD_LIBRARY(mitchell_shaeffer ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## NOBLE PURKINJE 1962 ##############################
set(MODEL_FILE_CPU "noble_1962.c")
set(COMMON_HEADERS model_common.h noble_1962.h)

ADD_LIBRARY(noble_1962 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## TEN TUSSHER 2006 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2006.c")
set(COMMON_HEADERS model_common.h ten_tusscher_2006.h)

ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## BEELER REUTER 1977 ##############################
set(MODEL_FILE_CPU "beeler_reuter_1977.c")
set(COMMON_HEADERS model_common.h beeler_reuter_1977.h)

ADD_LIBRARY(beeler_reuter_1977 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

