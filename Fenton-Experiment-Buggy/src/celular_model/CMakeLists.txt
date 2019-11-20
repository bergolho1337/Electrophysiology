set(SOURCE_FILES_CELL_MODEL celular_model.h celular_model.cpp)

ADD_LIBRARY(cell_model STATIC ${SOURCE_FILES_CELL_MODEL})

if(NOT MSVC)
    set_property(TARGET cell_model PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()