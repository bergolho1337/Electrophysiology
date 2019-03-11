set(SOURCE_FILES_MONODOMAIN monodomain.h monodomain.cpp)

ADD_LIBRARY(monodomain STATIC ${SOURCE_FILES_MONODOMAIN})

if(NOT MSVC)
    set_property(TARGET monodomain PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
