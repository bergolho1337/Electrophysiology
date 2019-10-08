cmake_minimum_required(VERSION 3.0.1)
project(FentonExperiment)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fno-strict-aliasing -Wall -Wno-unused-function -Wno-int-in-bool-context")

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    #message("${CMAKE_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)

#FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)

ADD_SUBDIRECTORY(src/monodomain)
ADD_SUBDIRECTORY(src/celular_model)
ADD_SUBDIRECTORY(src/config)
ADD_SUBDIRECTORY(src/plot)
ADD_SUBDIRECTORY(src/stimuli)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(FentonExperiment src/main.cpp)

ADD_DEPENDENCIES(FentonExperiment monodomain config plot stimuli cell_model utils)

#TARGET_LINK_LIBRARIES(FentonExperiment monodomain config plot stimuli cell_model utils dl m Eigen3::Eigen)
TARGET_LINK_LIBRARIES(FentonExperiment monodomain config plot stimuli cell_model utils dl m)

