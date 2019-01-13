cmake_minimum_required(VERSION 3.0.1)
project(stimulus)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

ADD_EXECUTABLE(stimulus stimulus.cpp)

TARGET_LINK_LIBRARIES(stimulus m)
