set(SOURCE_FILES_MONODOMAIN monodomain.h monodomain.cpp)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

ADD_LIBRARY(monodomain STATIC ${SOURCE_FILES_MONODOMAIN})

#if(NOT MSVC)
#    set_property(TARGET monodomain PROPERTY POSITION_INDEPENDENT_CODE ON)
#endif()
