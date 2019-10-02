SET(LIB_TYPE MODULE) #set the lib type

############## MITCHELL SHAEFFER 2003 ##############################
set(MODEL_FILE_CPU "mitchell_shaeffer_2003.c")
set(COMMON_HEADERS model_common.h mitchell_shaeffer_2003.h)

ADD_LIBRARY(mitchell_shaeffer_2003 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

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

############## DIFRANCESCO 1985 ##############################
set(MODEL_FILE_CPU "difrancesco_1985.c")
set(COMMON_HEADERS model_common.h difrancesco_1985.h)

ADD_LIBRARY(difrancesco_1985 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## ARPF 2009 ##############################
set(MODEL_FILE_CPU "stewart_aslanidi_noble_2009.c")
set(COMMON_HEADERS model_common.h stewart_aslanidi_noble_2009.h)

ADD_LIBRARY(stewart_aslanidi_noble_2009 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## LI & RUDY 2011 ##############################
set(MODEL_FILE_CPU "li_rudy_2011.c")
set(COMMON_HEADERS model_common.h li_rudy_2011.h)

ADD_LIBRARY(li_rudy_2011 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## TEN TUSSCHER ENDO 2004 ##############################
set(MODEL_FILE_CPU "ten_tusscher_endo_2004.c")
set(COMMON_HEADERS model_common.h ten_tusscher_endo_2004.h)

ADD_LIBRARY(ten_tusscher_endo_2004 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################

############## TEN TUSSCHER FIG 4b 2004 ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_Fig4b.c")
set(COMMON_HEADERS model_common.h ten_tusscher_3_Fig4b.h)

ADD_LIBRARY(ten_tusscher_3_Fig4b ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})

#########################################################
