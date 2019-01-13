cmake_minimum_required(VERSION 3.0.1)
project(CellModelSolver)

IF (CMAKE_VERSION VERSION_LESS "3.1")
    IF (CMAKE_C_COMPILER_ID STREQUAL "GNU")
        SET(CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    ENDIF(CMAKE_C_COMPILER_ID)
ELSE()
    SET(CMAKE_C_STANDARD 99)
ENDIF()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared_libs )

SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing -Wall -Wno-unused-function")

ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/solver)
ADD_SUBDIRECTORY(src/hash)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/string)

ADD_EXECUTABLE(CellModelSolver src/main.c)

ADD_DEPENDENCIES(CellModelSolver utils solvers hashes ini_parser string )

TARGET_LINK_LIBRARIES(CellModelSolver solvers utils hashes ini_parser string dl m)

ADD_SUBDIRECTORY(src/models_library)

SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.h)

ADD_SUBDIRECTORY(src/stimuli_library)
###

